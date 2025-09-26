***
*     Author:         Peter Dunn
*     Creation date:  15 February 1999
*     Last revision:  25 Septembeer 2025
*
* IN R, print statements for debugging as follows:
*     call dblepr("The value of y is ",-1, y, 1)
*     call intpr("The value of exact is ",-1, exact, 1)
******************************************************************
** Density-only routines
******************************************************************


******************************************************************
      subroutine twpdf(p, phi, y, mu, exact, verbose,
     &               funvalue, exitstatus, relerr, its )
*     Calculates the density of the log-likelihood
*     function of a Poisson-gamma distribution by inverting the
*     moment generating function.
*
*     IN:   p, phi, y, mu, exact, verbose
*     OUT:  funvalue, exitstatus, relerr, its
*     NOTE:  R shoudl always be supplying mu = 1

      double precision  p, phi, y, funvalue, mu, savemu,
     &                  lambda, calclambda, pi, area,
     &                  result, relerr, aimrerr
      integer  ier, maxit, m, iteratn, exitstatus,
     &         its, exact, verbose

* VARIABLES:
*    y        : the value at which the function is to be evaluated
*    x        : the range over which the integral is to be integrated;
*               an internal variable and NOT the value at which
*               the function is to be evaluated
*    c        : the constant mapping to the distribution with mean=1
*    lambda   : for 1<p<2, P(Y=0) = dexp( -lambda )
*    p        : the index (ie variance function is V(mu) = mu^p)
*    phi      : the dispersion parameter
*    funvalue : the value of the function at the given value of  x
*    bound    : The bound using Chebyshev theorm.
*    exitstatus:  1  if relative error is smaller than wished (aimrerr)
*                -1  if not, but the absolute error is less than aimrerr
*               -10  if neither rel or abs error any good
*    exact    : 1 if the exact zero acceleration algorithms is used;
*               0 if the approx zeros acceleration algorithm is used.
*    verbose  : 1 to print lots of diagnostic information
*             : 0 to keep quiet

*     Set defaults
      exitstatus = 1
      relerr = 0.0d00
      its = 0
      verbose = 0
      
*     SPECIAL CASES {
      if ( ( p .LT. 2.0d00 ) .AND. 
     &     ( p .GT. 1.0d00 ) ) then

         if ( y .LT. 0.0d00 ) then
            funvalue = 0.0d00

            return
         elseif ( y .EQ. 0.0d00 ) then
*           If  y = 0 in P-G case, then we don't need to integrate, just
*           evaluate; otherwise, integrate

            lambda = calclambda(p, phi, mu)
            funvalue = dexp( -1.0d00 * lambda )

            return
         endif
      elseif ( p .GE. 2.0d00 ) then

         if ( y .LE. 0.0d00 ) then
            funvalue = 0.0d00

            return
         endif
      endif

*     } END SPECIAL CASES
* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Specia cases handled in R. Remove for efficiency?

***

*     SET ACCURACY REQUIREMENTS
*     maximum number of iterations in calculating errors
      maxit = 100
      aimrerr = 1.0d-10

*     set other parameters
      m = -1
      pi = dacos( -1.0d00 )
      area = 0.0d00
      iteratn = 0
      relerr = 1.0d00
      ier = 0
      savemu = mu

*     For Tweedie distributions,
*        f(y; \mu, \phi) = cf(cy; c\mu, c^{2-p}\phi)
*     When mu\ne 1, we can set c = 1/\mu and then we
*     evaluate  f( y/\mu, 1, (\mu)^{2-p}\phi)
*     and the actual density should be c times this.
        
*     DETERMINE a(y,phi) (except in normal and gamma cases, where
*     the distribution is determined from the closed form result).

      if ( exact .EQ. 1 ) then
*        Use exact zero acceleration algorithm
         if ( (p .GT. 1.0d00).AND.(p .LE. 2.0d00) ) then
            call smallp( p, phi, y, mu, aimrerr, result,
     &                   maxit, ier, exitstatus, relerr, 
     &                   its, verbose )
*           This has found P(y | y>0); must convert back to P(Y)

         elseif ( p .GT. 2.0d00 ) then

            call bigp( p, phi, y, mu, aimrerr, result,
     &                 maxit, ier, exitstatus, relerr, 
     &                 its, verbose )

         endif

      else

*        Use approx zero acceleration algorithm
         call PDFapproxz( p, phi, y, mu, aimrerr, result,
     &                maxit, ier, exitstatus, relerr, 
     &                its, verbose )
*         result = -1000.0d00

      endif

      funvalue = result

*     Some tidy-ups
      if ( funvalue .LT. 0.0d00 ) funvalue = 0.0d00

      return
      end



******************************************************************
******************************************************************
      subroutine PDFapproxz( p, phi, y, mu, aimrerr, result,
     &              maxit, ier, exitstatus, relerr, its, verbose )
*     Calculates the density for distributions when  p > 2, using APPROX Sidi acceleration

      double precision  p, phi, y,  pi, area, aimrerr,
     &                  relerr, result, zero1, zero2,
     &                  mmatrix(2, 200), nmatrix(2, 200),
     &                  xvec(200), w, wold(3), area0, 
     &                  mu, f, f2, tmax, kmax, sumarea
      integer  its, ier, maxit, flag, exitstatus, itsidi, 
     &         mmax, verbose
      external  f, f2

*     VARIABLES:
*     exitstatus:    1  if relative error is smaller than wished (aimrerr)
*                   -1  if not, but the absolute error is less than aimrerr
*                  -10  if neither rel or abs error any good
*
*     SET OTHER PARAMETERS
      pi = dacos( -1.0d00 )
      area = 0.0d00
      area0 = 0.0d00
      its = 0
      relerr = 1.0d00
      flag = 0
      itsidi = 0

      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00

      if ( p .GT. 2.0d00 ) then

         kmax = 0.0d00
         tmax = 0.0d00
         mmax = -1

         if ( y .LT. 1.0d00 ) then
*           Find kmax
            call findkmax(p, phi, y, kmax, tmax, mmax, ier)
         endif

         zero2 = 0.0d00

*        Now, while t < tmax, do not use the Sidi acceleration. Instead, use ...

 400     if ( zero2 .LE. tmax ) then
*           get next zeros: jump the zeros by pi/(y)
            zero1 = zero2
            zero2 = zero2 + ( pi / y )

            call gaussq( f, result, zero1, zero2, 
     &                   p, phi, y, mu )
            area0 = area0 + result
            its = its + 1

            goto 400
         endif
         
         xvec(1) = zero2
      else
*        The initial region is done separately
*        since it isn't actually between zeros.
         zero1 = 0.0d00
         zero2 = pi / ( 2.0d00 * y )
         xvec(1) = zero2

         call gaussq( f2, area0, zero1, zero2, 
     &                p, phi, y, mu )
         its = its + 1

      endif

******Works OK to here

*     Now for Sidi acceleration!
*     Do some more integrations and use Sidi acceleration
  500 if (    ( ( itsidi .LT. 10 ) .AND.
     &          ( flag .NE. 1 )
     &        )
     &     .OR.
     &        ( ( itsidi .LT. maxit ) .AND.
     &          ( flag .NE. 1 ) .AND.
     &          ( abs(relerr) .GT. aimrerr )
     &        )
     &   ) then

*        get next zeros: jump the zeros by pi/(y)
         zero1 = zero2
         zero2 = zero2 + ( pi / y )

*        integrate between zeros
         if ( p .GT. 2.0d00 ) then
            call gaussq( f,  result, zero1, zero2, 
     &                   p, phi, y, mu )
         else
            call gaussq( f2, result, zero1, zero2, 
     &                   p, phi, y, mu )
         endif
         
*        Update interation count         
         its = its + 1
         itsidi = itsidi + 1

*        accelerate convergence of infinite sequence
         xvec( itsidi + 1 ) = zero2

         call accelerate( area, result, xvec, mmatrix,
     &        nmatrix, w, itsidi, relerr, wold, 
     &        sumarea, flag, verbose )

         relerr = ( abs( w - wold(1) ) + 
     &              abs( ( w - wold(2) ) ) ) / w 

         area = area + result

         goto 500

      endif

      result = ( area0 + w ) / pi
      if ( result .LT. 0.0d00 ) then
         result = 0.0d00
      endif
*     occasionally, a result may be very small, but negative

**     We have integrated to find  \int_0^{\infty}.  The integral
**     required actually goes from -infty to +infty, but is
**     symmetric about the y-axis, so the integral is _twice_ the result
**     obtained above.

      if ( flag .EQ. 1 ) then
         ier = -10
      endif

*     Determine the error
*     (Keep in this order so the most important aspect is returned)
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     we have good relative error, so that line should be OK.
      if ( abs(w) .LT. aimrerr ) then
         exitstatus = -1
      else
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) exitstatus = 1

      return
      end

****************************************************************
****************************************************************
      double precision function f(p, phi, y, mu, x)
*     A function to be numerically integrated

      double precision  x, p, phi, y, mu, rl, im

*   Tweedie EDM's are characterised by variance functions 
*      V( mu ) = phi * mu^p
*   for some p, where  p  is _outside_ the range (0, 1).
*
*   The function to be integrated is a Tweedie EDM of the form
*      f = b(y, phi) * exp{ (1/phi) * (y*(mu^[1-p]-1)/(1-p) -
*                                      (mu^[2-p])) /(2-p)}         (1)
*   which is of the form
*      f = a(y,phi) exp{ (1/phi) * (y * theta - kappa(theta) ) }  (2)
*
*   The density cannot be written in closed form when 1<p<2
*   due to the function outside the exponential
*   (ie, b(y, phi) in (1), or a(y,phi) in (2)),
*   and so there are two options:  (i) Evaluate this function by
*   an infinite summation; or (ii) evaluate by an infinitie integral, as
*   here.  This approach is prefered, since even when the density
*   can be written in closed form, the numerical evaluation
*   is still valid.
*
*   Given the cgf, K(.) for a distribution, the distribution can
*   be reconstructed:
*
*      f_Y(y) = 1/2*pi  \int_{\infty}^{\infty} dexp( K(ix) - ixy ) dx
*
*   This function then defines the integrand  dexp( K(ix) - ixy )
*   in the above expression.

* VARIABLES:
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   y          : the point at which the density function is to be evaluated
*   mu         : the mean value (here, a dummy 
*                to ensure gaussq  works correclt)
*   x          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf

      mu = 1
* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX remove mu for efficiency
*     mu  is listed an as unused input  

      call kt(p, phi, y, mu, x, rl, im)
      f = dexp( rl ) * dcos( im )

      return
      end

*******************************************************************
*******************************************************************
      double precision function f2( p, phi, y, mu, x )
*     A function to be numerically integrated in density, 1<p<2
*     using the conditional density.  Note that
*        lim (n->infty) dexp( Re k) = dexp( -lambda).
*
*     IN:  p, phi, y, mu, x
*     OUT: f2

* MAJOR VARIABLES:
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   y          : the point at which the density function is to be evaluated
*   mu         : the mean value
*   x          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
*

      double precision  x, p, phi, y, rl, im, mu,
     &                  calclambda, lambda

      mu = 1.0d00
      lambda = calclambda( p, phi, mu )

      if ( x .EQ. 0.0d00 ) then
         f2 = 1.0d00
*        True in the limit
      else
         call kt(p, phi, y, mu, x, rl, im)
         f2 = dexp( rl ) * dcos( im ) - 
     &         dexp( -lambda ) * dcos( x*y )
      endif

      return
      end

*****************************************************************
      subroutine smallp( p, phi, y, mu, aimrerr, result,
     &                   maxit, ier, exitstatus, relerr, its, 
     &                   verbose )
*     Calculates the density in the case 1 < p < 2
      double precision  p, phi, y,  pi, area, aimrerr,
     &         relerr, result, zero1, zero2,
     &         f, f2, g, intim, flo, fhi, t0, dk,
     &         mmatrix(2, 200), nmatrix(2, 200),
     &         xvec(500), w, wold(3), area0, 
     &         sumarea, sfzero2, mu, area1, 
     &         sfzero, zerofn, zerodfn, 
     &         lower, upper, tstep,     
     &         zarea0, z1lo, z1hi, zdelta, resultp
      integer  m, iteratn, ier, maxit, flag, numzr, tier,
     &         exitstatus, its, i, totalits,
     &         verbose, sbuffer, goflag

      external  f, g, intim, dk, sfzero2, zerofn,
     &          zerodfn, sfzero, f2

*     VARIABLES:
*     exitstatus:    1  if relative error is smaller than wished (aimrerr)
*                   -1  if not, but the absolute error is less than aimrerr
*                  -10  if neither rel or abs error any good

*     SET OTHER PARAMETERS
      m = -1
      pi = dacos( -1.0d00 )
      
      w = 0.0d00
      lambda = 0.0d00

      area = 0.0d00
      area0 = 0.0d00
      area1 = 0.0d00
      iteratn = 0
      relerr = 1.0d00
      sbuffer = 0
      flag = 0
      totalits = 1
*     totalits = 1 (not 0) as the first region is done separately; 
*     this is that one region     
      ier = 0
      tier = 0
      resultp = 0.0d00

      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00

*     FIND FIRST ZERO
      zero1 = 0.0d00
      
*     Find bounds on the other zero
      call findsp( p, mu, phi, y, lower, upper, flo, fhi )

*     This is linear interpolation between lower and upper:
      t0 = upper - fhi * ( upper - lower ) / ( fhi - flo )

      zero2 = sfzero( p, phi, y, lower, upper, t0, 
     &               zerofn, zerodfn, ier )
      tier = tier + ier

*     FIND FIRST AREA
*     The first region can be very strange.  For care, we use a
*     very high order Gaussian method, and break the region into
*     small pieces and operate on each separately.  Any funny business
*     should then hopefully be cornered.
      numzr = 20
* WAS 20:  Changed 07/Dec/2005
      
      if ( verbose .EQ. 1 ) then
         call dblepr(" Integrating between ",-1,zero1,1)
         call dblepr(" and ",-1,zero2,1)
         call intpr(" using this many sub-regions: ",-1,numzr,1)
      endif
      
      
      zdelta = zero2 / dble( numzr )
      z1lo = 0.0d00
      z1hi = 0.0d00
      do i = 1, numzr
         zarea0 = 0.0d00
         z1lo = z1hi
         z1hi = z1hi + zdelta
         call gaussq(f2, zarea0, z1lo, z1hi, p, phi, y, mu )
         area0 = area0 + zarea0
      enddo
*     So that is one iteration (between t=0 and t=<first zero>: 
*     hence totalits = 1 up to here

      zero1 = zero2
      tstep = zero2 / 2.0d00

*     NOW, DO A FEW MORE ITS FOR SAFETY'S SAKE
*     as the regions can be dodgy when 1<p<2
      goflag = 1
      flag = 0
      area1 = 0.0d00

550   if ( goflag .EQ. 1 ) then

         totalits = totalits + 1

*        FIND THE ZERO
*        We don't jump too far ahead of ourselves--especially
*        early on, when some intervals can be pretty dodgy.
         lower = zero1 + tstep * 0.05d00
         upper = zero1 + 0.3d00 * tstep

         flo = zerofn( p, phi, y, lower )
         fhi = zerofn( p, phi, y, upper )

*        Try harder to bound the zero        
 650     if ( ( flo * fhi ) .GT. 0.0d00 ) then

            lower = upper
            upper = upper + 0.5d00*tstep

            flo = zerofn(p, phi, y, lower)
            fhi = zerofn(p, phi, y, upper)

            goto 650
         
         endif


         zero2 = sfzero( p, phi, y, lower, upper, t0,
     &                   zerofn, zerodfn, ier )
         tier = tier + ier

*        Keep track of last result, too
         resultp = result

*        INTEGRATE
         call gaussq( f2, result, zero1, zero2, 
     &                p, phi, y, mu )
         if ( verbose .EQ. 1 ) then
            call dblepr(" Integrating between ",-1,zero1,1)
            call dblepr(" and ",-1,zero2,1)
         endif

*        SUM AREA
         area1 = area1 + result

*        PREPARE FOR NEXT ITERATION
         tstep = zero2 - zero1
         zero1 = zero2 
         t0 = zero2 + ( 0.8d00 * tstep )

*        See if we can stop now
         if ( sbuffer .GE. 3 ) then
            goflag = 0
         endif

         sbuffer = sbuffer + 1

*        Also stop if the areas are getting really small
*         if ( ( abs(result) .LT. aimrerr/10000.d00 ) .AND.
*     &        (abs(resultp) .LT. aimrerr/10000.d00 ) .AND.
*     &        ( its .GT. 10 ) ) then
*            goflag = 0
*         endif

*         if ( result .EQ. 0.0d00 ) then
*            go=0
*         endif
* Those lines above commented out 15 Sep 2005

         goto 550

      endif

*     NOW, INTEGRATE WITH ACCELERATION
*     We only need to do this if the regions have any significant area
*      if ( ( abs(result) .LT. aimrerr/10000.d00 ) .AND.
*     &     (abs(resultp) .LT. aimrerr/10000.d00 ) .AND.
*     &     ( its .GT. 10 ) ) then
*         goflag = 0
*         exitflag = 1
*      else
*         goflag = 1
*      endif
      goflag = 1
      flag = 0
      its = 0
      area = 0.0d00
      xvec(1) = zero2

 1550 if ( goflag .EQ. 1 ) then

         its = its + 1
         totalits = totalits + 1
         
*        FIND THE ZERO
*        Our jumps here can be a little more bold, since most of the
*        initial antics should have been sorted out.
         lower = zero1 + 0.05d00 * tstep
         upper = zero1 + 0.8d00 *tstep
 
*        Note:  zero1 has been set above

         flo = zerofn( p, phi, y, lower )
         fhi = zerofn( p, phi, y, upper)

 1650    if ( ( flo * fhi ) .GT. 0.0d00 ) then

            lower = upper
            upper = upper + 0.5d00*tstep

            flo = zerofn(p, phi, y, lower)
            fhi = zerofn(p, phi, y, upper)

            goto 1650

         endif

*        This is linear interpolation:
         t0 = lower - flo * (upper -  lower ) / 
     &            ( fhi - flo )  
         zero2 = sfzero( p, phi, y, lower, upper, t0,
     &                   zerofn, zerodfn, ier )
          tier = tier + ier

            call gaussq( f2, result, zero1, zero2,
     &                     p, phi, y, mu )

*           ACCELERATE CONVERGENCE
            xvec( its + 1 ) = zero2
            call accelerate( area, result, xvec, mmatrix,
     &           nmatrix, w, its, relerr, wold, sumarea, 
     &           flag, verbose )
            
            if ( its .GE. 3 ) then
               relerr =  ( abs(   w - wold(1) ) + 
     &                     abs( ( w - wold(2) ) ) )
     &               / (area0 + area1 + w)
            endif
         
         if ( flag .EQ.1 ) then
*            print *,'Machine limits being reached...'
         endif

*        SUM AREA
         area = area + result

*        PREPARE FOR NEXT ITERATION
         tstep = zero2 - zero1
         zero1 = zero2 
         t0 = zero2 + ( 0.8d00 * tstep )

*        NOTE IF FLAG=1, the limits of the machine are being reached.
*        We stop and report that required accuracy may not be achieved.

         if (     ( its .LT. 3 ) 
     &         .OR.
     &            ( ( its .LT. maxit ) .AND.
     &              ( abs(relerr) .GT. aimrerr )
     &            )
     &      ) then

*           THEN keep going...
            goflag = 1
         else
*           THEN stop
            goflag = 0
         endif

*        Now sometimes we get a very small w, whose relative
*        error isn't that small.  Since w is small, we can make
*        do with it anyway.  So we check this case too.

         if ( flag .EQ. 1 ) then
            ier = -90
            tier = tier + ier
         endif

         goto 1550
      endif

      result = area0 + area1 + w
      result = result / pi
      ier = tier
      
*     Now report the total number of iterations.
*     This is stored as  totalits; we used  its  above
*     as it was needed for Sidi acceleration.  Now we
*     revert to  its
      its = totalits

*     Recall that the function we integrate (which one?) has the division by
*     (1.0d00 - dexp( -lambda ) ) already built into it.
            
*     Determine the error
*     In this order, the most important aspect is returned.
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     we have good relative error, so that line should be OK.
      if ( ( abs(w - wold(1)) + abs(w - wold(2)) ) .LT. aimrerr ) then
*        Absolute error isn't too bad
         exitstatus = -1
      else
*        Even absolute error isn't too good
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) then
*        Relative error is inside the required accuracy
         exitstatus = 1
      endif

      return

      end

*****************************************************************
****************************************************************
      subroutine bigp( p, phi, y, mu, aimrerr, result,
     &                    maxit, ier, exitstatus, relerr, 
     &                    its, verbose )
*     Calculates the density in the case p > 2.

      double precision  p, phi, y, mu, pi, area, aimrerr,
     &          relerr, result, zero1, zero2, zero,
     &          f, g, zlo, zhi, intim, flo, fhi, kmax,
     &          tmax, mmatrix(2, 200), nmatrix(2, 200),
     &          xvec(200), w, wold(3), area0, sumarea,
     &          dk, sfzero2, diff, largest, smallest
      integer  m, its, mmax, firstm, ier, maxit, flag,
     &         allok, kmaxok, tier, exitstatus, verbose
      external  f, g, intim, dk, sfzero2

*     VARIABLES:
*     exitstatus:    1  if relative error is smaller than wished (aimrerr)
*                   -1  if not, but the absolute error is less than aimrerr
*                  -10  if neither rel or abs error any good

      if ( verbose .EQ. 1 ) then
         call dblepr("Using p>2 code since p = ",-1, p, 1)
      endif

*     SET ACCURACY REQUIREMENTS
      largest = 1.0d30
      smallest = 1.0d-30

*     SET OTHER PARAMETERS
      m = -1
      pi = dacos( -1.0d00 )
      area = 0.0d00
      area0 = 0.0d00
      its = 0
      relerr = 1.0d00
      flag = 0
      tier = 0
      allok = 1

      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00

      if ( y .GE. 1.0d00 ) then
      
         if ( verbose .EQ. 1 ) then
            call dblepr(" Using y .GE. 1 since y = ",-1, y, 1)
         endif
*        In this case, Im(k) heads down straight away.

*        FIND ZEROS
         zero1 = 0.0d00

*        An approximation to the first zero:
         zero = pi / ( 2.0d00 * y )

*        Bracket first zero
         zlo = 0.9d00 * pi / (2.0d00 * y )
         
         if ( y .GT. 1.0d00 ) then
            zhi = pi / (2.0d00 * ( y - 1.0d00 ) )
            fhi = intim( p, phi, y, zhi, m )
         else
            zhi = zero * 2.0d00  
            fhi = intim( p, phi, y, zhi, m )
         endif
         flo = intim( p, phi, y, zlo, m )
         
         allok = 1
 
 565     if ( ( allok .EQ. 1 ) .AND.
     &        (fhi * flo ) .GT. 0.0d00 ) then

            zlo = zhi 
            zhi = zhi * 1.5d00

            flo = intim( p, phi, y, zlo, m )
            fhi = intim( p, phi, y, zhi, m )

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            goto 565

         endif

         if ( zhi .GT. largest/10.0d00 ) allok = 0

         if ( allok .EQ. 0 ) then
            ier = -50
            tier = tier + ier
            result = 0.0d00
            exitstatus = -10
            return
         endif

         if ( verbose .EQ. 1 ) then
            call dblepr("  First zero bound by ",-1, zlo, 1)
            call dblepr("  and ",-1, zhi, 1)
         endif
         zero2 = sfzero2( p, phi, y, zlo, zhi, zero, 
     &                   intim, dk, m, ier )
         if ( verbose .EQ. 1 ) then
            call dblepr("  Found first zero = ",-1, zero2, 1)
         endif

         xvec( 1 ) = zero2

*        special case: between 0 and first zero
         if ( verbose .EQ. 1 ) then
            call dblepr("   Integrating between y = ",-1, zero1, 1)
            call dblepr("                   and y = ",-1, zero2, 1)
         endif
         call gaussq( f, area0, zero1, zero2, 
     &                p, phi, y, mu )
         if ( verbose .EQ. 1 ) then
            call dblepr("   Area is ",-1, area0, 1)
         endif

*        NOTE:  We don't update the iteration count here, since
*        we keep it all to work with Sidi acceleration; we instead
*        add one later (after accelerating)

*        Now do some more integrations and use Sidi acceleration
!   500    if (    ( ( its .LT. 4 ) .AND.
!      &             ( flag .NE. 1 )
!      &           )
!      &        .OR.
!      &           ( ( its .LT. maxit ) .AND.
!      &             ( flag .NE. 1 ) .AND.
!      &             ( abs(relerr) .GT. aimrerr )
!      &           ) 
!      &      ) then
  500    if (    ( its .LT. 4 )
     &        .OR.
     &           ( ( its .LT. maxit ) .AND.
     &             ( abs(relerr) .GT. aimrerr )
     &           ) 
     &      ) then
     
            if ( verbose .EQ. 1 ) then
               call intpr("   Iterating; iteration ",-1, its, 1)
            endif
     

*           get next zeros
            m = m - 1
            zero1 = zero2
            zero = zero2
            zlo = zero2
            zhi = zero2 * 1.5d00

            if ( zhi .GT. largest/10.0d00 ) then

               allok = 0
               flo = 0.0d00
               fhi = 0.0d00

            else

               allok = 1
               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )
 
            endif

 765        if ( ( allok .EQ. 1 ) .AND.
     &           (fhi * flo ) .GT. 0.0d00 ) then

               zlo = zhi
               zhi = zhi * 1.5d00

               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

*              Linearly interpolate for zero:
               zero = zhi - fhi * ( zhi - zlo ) / 
     &                  ( fhi - flo )

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               goto 765

            endif

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10
               return
            endif
     
            zero2 = sfzero2( p, phi, y, zlo, zhi, zero, 
     &                      intim, dk, m, ier )

            if ( ier .NE. 0 ) then
               tier = tier + ier
            endif


*           integrate between zeros
            if ( verbose .EQ. 1 ) then
               call dblepr("   Integrating between y = ",-1, zero1, 1)
               call dblepr("                   and y = ",-1, zero2, 1)
            endif

            call gaussq( f, result, zero1, zero2, 
     &                   p, phi, y, mu )
            if ( verbose .EQ. 1 ) then
               call dblepr("   Giving area = ",-1, result, 1)
            endif

*           Update iteration count
            its = its + 1

*           accelerate convergence of infinite sequence
            xvec( its+1 ) = zero2

            call accelerate( area, result, xvec, mmatrix,
     &           nmatrix, w, its, relerr, wold, sumarea, 
     &           flag, verbose )
            if ( verbose .EQ. 1 ) then
               call dblepr("   Accelerating; w = ",-1, w, 1)
            endif

            relerr = (abs( w - wold(1)) + abs(( w - wold(2))))
     &               / (area0 + w)       

            area = area + result

            if ( verbose .EQ. 1 ) then
               call dblepr("  Area = ",-1, area, 1)
            endif
            
            goto 500

         endif

         if ( ( its .GE.  maxit ) .AND.
     &        ( abs(relerr) .GT. aimrerr) ) then
            ier = -40
            tier = tier + ier
         endif
         if ( flag .EQ. 1 ) then
            ier = -70
            tier = tier + ier
         endif

         result = area0 + w

*        Now, the very first integration has not been counted
*        (since that can foul up the Sidi acceleration iteration count)
*        so update now
         its = its + 1

      else
*        that is:    if ( y .LT. 1) then
*        In this case, Im(k) may head up before going to  -infinity
         if ( verbose .EQ. 1 ) then
            call dblepr(" Using y .LT. 1 since y = ",-1, y, 1)
         endif


*        FIND k_max AND t_max
         kmaxok = 1
         call findkmax(p, phi, y, kmax, tmax, mmax, ier)
         if ( verbose .EQ. 1 ) then
            call dblepr(" Found kmax = ",-1, kmax, 1)
         endif

         if ( ier .NE. 0 ) then
            tier = tier + ier
            kmaxok = 0
         endif
         

         if ( kmax .LT. pi/2.0d00 ) then

            if ( verbose .EQ. 1 ) then
               call dblepr(" Using kmax .LT. pi/2; kmax = ",
     &                     -1, kmax, 1)
            endif
            m = -1

*           FIND ZEROS
            zero1 = 0.0d00
            zero = tmax + pi/( 2.0d00*y )

*           BOUNDS ON `OTHER' ZERO:
            zlo = tmax
            zhi = zero*2.0d00

            if ( zhi .GT. largest/10.0d00 ) then

               allok = 0
               flo = 0.0d00
               fhi = 0.0d00

            else

               allok = 1
               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

            endif

 1565       if ( (allok .EQ. 1 ) .AND.
     &            ( (fhi * flo ) .GT. 0.0d00 ) ) then

               zlo = zhi
               zhi = zhi * 1.5d00

               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               goto 1565

            endif

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10
               return
            endif

            zero2 = sfzero2( p, phi, y, zlo, zhi, zero, 
     &                      intim, dk, m, ier )

            if ( ier .NE. 0 ) then
               tier = tier + ier
            endif

            xvec(1) = zero2

*           integrate between zeros
            if ( verbose .EQ. 1 ) then
               call dblepr("  Integrating between y = ",
     &                     -1, zero1, 1)
               call dblepr("                  and y = ",
     &                     -1, zero2, 1)
            endif
            
            call gaussq( f, area0, zero1, zero2, 
     &                   p, phi, y, mu )
            
            if ( verbose .EQ. 1 ) then
               call dblepr("  Giving area = ",-1, area0, 1)
            endif

*        NOTE:  We don't update the iteration count here, since
*        we keep it all to work with Sidi acceleration; we instead
*        add one later (after accelerating)


  600       if (    ( its .LT. 4 )
     &            .OR.
     &               ( ( its .LT. maxit ) .AND.
     &                 ( abs(relerr) .GT. aimrerr )
     &               ) ) then

            if ( verbose .EQ. 1 ) then
               call intpr("  Iteration ",-1, its, 1)
            endif

*              get next zeros
               m = m - 1
               diff = zero2 - zero1
               zero1 = zero2

               zlo = zero2 - 0.01d00 * diff
               zhi = zero2 + 2.0d00 * diff

               if ( zhi .GT. largest/10.d00 ) then

                  allok = 0
                  flo = 0.0d00
                  fhi = 0.0d00

               else

                  allok = 1
                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

               endif

 1665          if ( ( allok .EQ. 1 ) .AND.
     &              (fhi * flo ) .GT. 0.0d00 ) then

                  zlo = zhi
                  zhi = zhi * 1.5d00

                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

                  if ( zhi .GT. largest/10.0d000 ) then
                     allok = 0
                  endif

                  goto 1665
   
               endif

               if ( zhi .GT. largest/10.d00 ) allok = 0

               if ( allok .EQ. 0 ) then
                  ier = -50
               tier = tier + ier
                  result = 0.0d00
               exitstatus = -10
   
                  return

               endif

*              Approximate zero with linear interpolation
               zero = zlo - flo * ( zhi - zlo ) / 
     &                     ( fhi - flo )

               zero2 = sfzero2( p, phi, y, zlo, zhi, zero,
     &                         intim, dk, m, ier )

               if ( ier .NE. 0 ) then
                  tier = tier + ier
               endif


*              integrate between zeros
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Integrating between y = ",-1, zero1, 1)
                  call dblepr("                  and y = ",-1, zero2, 1)
               endif
               call gaussq( f, result, zero1, zero2, 
     &                      p, phi, y, mu )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Giving area ",-1, result, 1)
               endif
     
*              Update interation count
               its = its + 1

*              accelerate convergence of infinite sequence
               xvec( its+1 ) = zero2

               call accelerate( area, result, xvec, mmatrix, 
     &              nmatrix, w, its, relerr, wold, 
     &              sumarea, flag, verbose )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Accelerating; w = ",-1, w, 1)
               endif
               relerr = ( abs( w - wold(1) ) + 
     &                    abs( ( w - wold(2) ) ) )
     &                  / (area0 + w)       

               area = area + result
               
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Area = ",-1, area, 1)
               endif

               goto 600

            endif

            if ( ( its .GE.  maxit ) .AND.
     &           ( abs(relerr) .GT. aimrerr) ) then
               ier = -40
               tier = tier + ier
            endif

            if ( flag .EQ. 1 ) then
               ier = -70
               tier = tier + ier
            endif
      
            result = area0 + w

*           Now, the very first integration has not been counted
*           (since that can foul up the Sidi acceleration iteration count)
*           so update now
            its = its + 1

         else
         
*           that is:   case where kmax >= pi/2
*           In this case, the upward trend is goes above  pi/2, and so
*           the first zero will be at  m=0.

            if ( verbose .EQ. 1 ) then
               call dblepr(" Using kmax .GE. pi/2; kmax = ",
     &                     -1, kmax, 1)
            endif

*           Now, kmax may not have been found accurately.  IF, however,
*           the corresponding max > maxit, it won't matter a great deal
*           and we can proceed.  If not, accuracy cannot be ensured
*           unless the maximum  m  used is less than maxit.

            if ( ier .EQ. -80 ) then
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10
               return
            endif

            wold(1) = 0.0d00
            wold(2) = 0.0d00
            wold(3) = 0.0d00

*           FIND ZEROS
            zero1 = 0.0d00
            zero = pi / ( 2.0d00*( 1.0d00 - y ) )

            m = 0
            firstm = 1

            zlo = smallest
            zhi = tmax

            if ( zhi .GT. largest/10.d00 ) then
              allok = 0
              flo = 0.0d00
              fhi = 0.0d00
            else
               allok = 1
               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )
            endif

            diff = zhi - zlo

 2565       if ( ( allok .EQ. 1 ) .AND. 
     &           (fhi * flo ) .GT. 0.0d00 ) then

               zlo = zhi
               zhi = zhi + 0.1d00*diff

               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

               if ( zhi .GT. largest / 10.0d00 ) then
                  allok = 0
               endif

               goto 2565

            endif

            if ( zhi .GT. largest/10.d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10

               return
            endif

            zero2 = sfzero2( p, phi, y, zlo, zhi, zero, 
     &                      intim, dk, m, ier )

            if ( ier .NE. 0 ) then
               tier = tier + ier
            endif

            xvec( 1 ) = zero2

            if ( verbose .EQ. 1 ) then
               call dblepr("  Integrating between y = ",-1, zero1, 1)
               call dblepr("                  and y = ",-1, zero2, 1)
            endif
            call gaussq( f, area0, zero1, zero2, 
     &                   p, phi, y, mu )
            if ( verbose .EQ. 1 ) then
               call dblepr("  Giving area = ",-1, area0, 1)
            endif
*           NOTE:  We don't update the iteration count here, since
*           we keep it all to work with Sidi acceleration; we instead
*           add one later (after accelerating)


            diff = zero2 - zero1

  700       if (    ( its .LT. 4 )
     &           .OR.
     &              ( ( its .LT. maxit ) .AND.
     &                ( abs(relerr) .GT. aimrerr )
     &              ) ) then

               if ( verbose .EQ. 1 ) then
                  call intpr("  Iteration ", 1, its, 1)
               endif

*              get next zeros

               zlo = zero2 - 1.0d-05*diff
               zhi = zero2 + 2.0d00*diff

               zero1 = zero2

*              FIND THE NEXT VALUE OF m
               call nextm( tmax, mmax, zero2, m, firstm, 
     &                     zlo, zhi, zero )

               if ( zhi .GT. largest/10.d00 ) then

                  allok = 0
                  flo = 0.0d00
                  fhi = 0.0d00

               else

                  allok = 1
                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

               endif

 2665          if ( ( allok .EQ. 1 ) .AND.
     &              (fhi * flo ) .GT. 0.0d00 ) then

                  zlo = zhi
                  zhi = zhi * 1.5d00

                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

                  if ( zhi .GT. largest/10.0d00 ) then
                     allok = 0
                  endif

                  goto 2665

               endif

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               if ( allok .EQ. 0 ) then
                  ier = -50
               tier = tier + ier
                  result = 0.0d00
                  exitstatus = -10
   
                  return
               endif

               zero2 = sfzero2( p, phi, y, zlo, zhi, zero,
     &                         intim, dk, m, ier )

               if ( ier .NE. 0 ) then
                 tier = tier + ier
               endif

               if ( verbose .EQ. 1 ) then
                  call dblepr("  Iteragrating between y = ", 
     &                        1, zero1, 1)
                  call dblepr("                   and y = ", 
     &                        1, zero2, 1)
               endif
               call gaussq( f, result,  zero1, zero2,
     &                      p, phi, y, mu )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  giving area = ", 1, 
     &                        result, 1)
               endif

*              Update iteration count
               its = its + 1
               
*              accelerate convergence of infinite sequence
               xvec( its+1 ) = zero2
               call accelerate( area, result, xvec, mmatrix, 
     &            nmatrix, w, its, relerr, wold, sumarea, 
     &            flag, verbose )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Accelerating; w = ", 
     &                        1, w, 1)
               endif

               area = area + result
               
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Area = ", 1, area, 1)
               endif

               goto 700

            endif

            if ( ( its .GE.  maxit ) .AND.
     &           ( abs(relerr) .GT. aimrerr) ) then
               ier = -40
               tier = tier + ier
            endif
            if ( flag .EQ. 1 ) then
               ier = -70
               tier = tier + ier
            endif

            if ( kmaxok .EQ. 0 ) then
*              IF finding kmax didn't converge...

               if ( m .LT. mmax-1 )  kmaxok = 1
*              All should be OK if greatest value of m used
*              is less than the `turning' m value.
            endif

            if ( kmaxok .EQ. 0 )  then
               ier = -60
               tier = tier + ier
            endif

            result = area0 + w
*           The first integration has not been counted (that 
*           can foul up the Sidi acceleration) so update now
            its = its + 1
         endif
      endif

**     We have integrated to find  \int_0^{\infty}.  The integral
**     required actually goes from -infty to +infty, but is
**     symmetric about the y-axis, so the integral is _twice_ the result
**     obtained above.

      result = abs( result / pi )
*     occasionally, a result may be very small, but negative

      if ( flag .EQ. 1 ) then
         ier = -10
      endif

      if ( ier .NE. 0 ) then
               tier = tier + ier
      endif

*     Determine the error
*     (Keep in this order so the most important aspect is returned)
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     We have good relative error, so that line should be OK.
      if ( abs(w) .LT. aimrerr ) then
         exitstatus = -1
      else
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) exitstatus = 1

      if ( verbose .EQ. 1 ) then
         call dblepr("Final result = ", 1, result, 1)
      endif
      
      return
      end

*****************************************************************
*******************************************************************

      double precision function sfzero( p, phi, y, x1, x2, x0,
     &                                 fun, dfun, ier )
*     Uses modified Newton's Method to find a root between
*     x1 and x2 to find kmax

      double precision  df, dx, dxold, f, fh, fl, temp, 
     &                  xh, xl, x1, x2, fun, dfun, y, p, 
     &                  phi, x0
      integer  j, maxit, ier
      external  fun, dfun

*     SET PARAMETERS
      ier = 0
      maxit = 100
      
      fl = fun(p, phi, y, x1)
      fh = fun(p, phi, y, x2)

      if ( fl .EQ. 0.0d00 ) then
         sfzero = x1
         return
      elseif ( fh .EQ. 0.0d00 ) then
         sfzero = x2
         return
      elseif ( fl .LT. 0.0d00 ) then
         xl = x1
         xh = x2
      else
         xl = x2
         xh = x1
      endif
         
      sfzero = x0
      dxold = abs( x2-x1 )
      dx = dxold

      f = fun( p, phi, y, sfzero )
      df = dfun( p, phi, y, sfzero )

      do j = 1, maxit

         if (( (sfzero - xh) * df - f) * ((sfzero - xl) * df - f) 
     &                                              .GT. 0.0d00
     &                .OR.
     &        abs(2.0d00 * f) .GT. abs( dxold * df)) then

            dxold = dx
            dx = 0.5d00 * ( xh - xl )
            sfzero = xl + dx      
            if ( xl .EQ. sfzero ) return

         else

            dxold = dx
            if ( df .EQ. 0.0d00 ) return
               
            dx = f/df
            temp = sfzero
            sfzero = sfzero - dx
            if ( temp .EQ. sfzero ) return

         endif

!           if ( abs( dx ) .LT. 1.0d-13 ) return

         f = fun( p, phi, y, sfzero )
         df = dfun( p, phi, y, sfzero )

         if ( f .LT. 0.0d00 ) then
            xl = sfzero
         else
            xh = sfzero
         endif

      enddo

      ier = -20

      return
      end
      
      
*******************************************************************
*******************************************************************
      double precision function imgcgf(p, phi, y, x)
*     Calculates the imag part of the cgf

      double precision  x, p, phi, y, rl

      mu = 1.0d00
      call kt( p, phi, y, mu, x, rl, imgcgf )

      return
      end

*******************************************************************
*******************************************************************
      double precision function sfzero2( p, phi, y, x1, x2,
     &                       x0, fun, dfun, m, ier )
*     Uses a modified Newton's Method to find a root between
*     x1 and x2 to find kmax

      double precision  df, dx, dxold, f, fh, fl, temp, 
     &         xh, xl, x1, x2, fun, dfun, y,
     &         p, phi, x0
      integer  j, maxit, m, ier
      external  fun, dfun

*     SET PARAMETERS
      maxit = 100
      ier = 0

*         call dblepr("x1:",-1, x1, 1)
*         call dblepr("x2:",-1, x2, 1)
         
      fl = fun(p, phi, y, x1, m)
      fh = fun(p, phi, y, x2, m)

*      if ( ( fl .GT. 0.0d00 .AND. fh .GT. 0.0d00 ) 
*     &          .OR.
*     &     ( fl .LT. 0.0d00 .AND. fh .LT. 0.0d00 ) ) then
*
*         print *,'Error: Root must be bounded.'
*         print *,'xs are : ',x1,x2
*         print *,'giving :', fl, fh
*         print *,'for m = ',m
*         return
*
*      endif

      if ( fl .EQ. 0.0d00 ) then
      
         sfzero2 = x1
*      print *,'POI Return: sfzero2=',sfzero2
         return
      elseif ( fh .EQ. 0.0d00 ) then
         sfzero2 = x2
*      print *,'ASP Return: sfzero2=',sfzero2
         return
      elseif ( fl .LT. 0.0d00 ) then
         xl = x1
         xh = x2
      else
         xl = x2
         xh = x1
      endif
         
      if ( ( x0 .GT. xl ) .AND. ( x0 .LT. xh ) ) then
         sfzero2 = x0
      else
         sfzero2 = ( xl + xh ) / 2.0d00
      endif
      dxold = abs( x2 - x1 )
      dx = dxold

      f = fun( p, phi, y, sfzero2, m )
      df = dfun( p, phi, y, sfzero2 )

      do j = 1, maxit

         if ( ( ( sfzero2 - xh * df - f) * (sfzero2 - xl) * df - f) 
     &          .GT. 0.0d00
     &                .OR.
     &        abs(2.0d00 * f) .GT. abs( dxold * df) ) then
*           Then use bisection

            dxold = dx
            dx = 0.5d00 * ( xh - xl )
            sfzero2 = xl + dx
            if ( xl .EQ. sfzero2 ) then
               return
            endif
         else
*           Then use Newton's method
            dxold = dx
            dx = f/df
            temp = sfzero2
            sfzero2 = sfzero2 - dx
            if ( temp .EQ. sfzero2 ) then
               return
            endif
         endif
         
*         call dblepr("sfzero2:",-1, sfzero2, 1)
*         call intpr("Iteration:",-1, j, 1)

         if ( abs( dx ) .LT. 1.0d-11 ) then
            return
         endif

         f = fun( p, phi, y, sfzero2, m )
         df = dfun( p, phi, y, sfzero2 )

         if ( f .LT. 0.0d00 ) then
            xl = sfzero2
         else
            xh = sfzero2
         endif

      enddo

      ier = -30

      return
      end


*******************************************************************
*****************************************************************
      double precision function othzero(p, phi, y )
*     Finds the `other zero' in the density integrand
*     in finding a starting point.  (This `other root' is the zero
*     of dsin( Im k).)

*     IN:  p, phi, y, mu
*     OUT: othzero

      double precision  p, phi, y, pi, psi, inflec, tlo,
     &                  thi, t0, intim, dk, tmax, kmax,
     &                  sfzero2, smallest, largest,
     &                  flo, fhi, zstep
      integer  ier, mmax, maxit, m
      external  intim, dk, sfzero2

      largest = 1.0d30
      smallest= 1.0d-30

      ier = 0
      maxit = 100
      pi = dacos( -1.0d00 )

*     This is the point of inflection: a starting point for the zeroing:
      psi = ( pi / 2.0d00 ) * ( 1.0d00 - p ) / 
     &          ( 2.0d00 * p - 1.0d00 )
      inflec = datan( psi ) / ( ( 1.0d00 - p) * phi )

*     First, we establish whether Im k is going up or down
*     initially.  This established what value of  m  we are
*     using to get the zero.  If y>mu, we are heading up first.

      if ( y .GE. 1.0d00 ) then
*        In this case, we initially head down.
         m = -1
         kmax = 0.0d00
         tmax = 0.0d00
         tlo = 1.0d-05
         thi = inflec
      else
*        In this case, we have to find k_max.  To do so, we
*        solve Im k'(t) = 0.
         call findkmax( p, phi, y, 
     &                 kmax, tmax, mmax, ier )

         if ( kmax .GE. pi/2.0d00 ) then
            m = 0
            tlo = smallest
            thi = tmax
         else
            m = -1 
            tlo = min(tmax, inflec)
            thi = max(tmax, inflec)
         endif
      endif

      flo = intim(p, phi, y, tlo, m)
      fhi = intim(p, phi, y, thi, m)
      zstep = thi - tlo

 565  if ( ( flo * fhi ) .GT. 0.0d00 ) then
           tlo = thi
           thi = thi + 0.2d00 * zstep

           flo = intim(p, phi, y, tlo, m)
           fhi = intim(p, phi, y, thi, m)

           goto 565
        endif

*     Having established  m, this means we can find the root
*     we are after a bit more easily.

*     Linear interpolation for the estimate:
      t0 = tlo - flo * ( thi - tlo ) / ( fhi - flo )

      othzero =  sfzero2( p, phi, y, tlo, thi, t0,
     &                intim, dk, m, ier)

      return
      end

*****************************************************************
****************************************************************
      subroutine findkmax(p, phi, y, kmax, tmax, mmax, ier)
*     Finds  k_max and t_max; solve k'(t) = 0 with Newton's method.

      double precision  kmax, tmax, imgdcgf, kddasht, 
     &            imgcgf, p, largest, dpmmax, zhi, zlo, 
     &            phi, y, pi, psi, front, inner, z,flo, 
     &            fhi, dz1, dz2, z1, z2, sfzero, dk, dz
      integer  mmax, ier, allok, lrgint
      external  dk, kddasht, imgdcgf

      pi = dacos( -1.0d00 )
      ier = 0
      allok = 1
      largest = 1.0d30
      lrgint = 100000000

*     A good starting point is sometime crucial, so we spend a 
*     little time finding a decent one.

*     First, the point of inflection, which will always
*     work, but may be very slow.

      psi = ( pi / 2.0d00 ) * ( 1.0d00 - p ) / 
     &      ( 2.0d00*p - 1.0d00 )

*     psi is point of inflection, so start Newton's method at this point
      z2 = 1.0d00 / ( phi * (1.0d00 - p) ) * dtan( psi )
      
      dz2 = imgdcgf(p, phi, z2) - y

*     But an alternative, closer to the real value for small y 
*     can be found provided p > 2:
      if ( p .GT. 2.0d00) then
         front = -1.0d00 / ( phi * (1.0d00 - p ) )
         inner = ( 1.0d00 / y ) * 
     &             dcos( -pi / (2.0d00 * (1.0d00 - p)) )
*        NOTE:   inner  will always be positive when p>2

         z1 = front * ( inner ) ** ( p - 1.0d00 )
         dz1 = imgdcgf(p, phi, z1) - y

*     Now choose the starting point that is the larger, but is still
*     to the left of the k_max; this can be done by examining the
*     derivative at each of the two points.  
*     In effect, we choose the largest starting value for which
*     the derivative is positive.

         if ( dz1 .GT. 0.0d00 ) then
*           If both starting points are to the left of kmax,
*           choose the larger
            if ( z1 .GT. z2 ) then
               z = z1
               dz = dz1
            else
               z = z2
               dz = dz2
            endif
         else
            if ( dz2 .LT. 0.0d00 ) then
*              That is, both to the right; so, use the smaller
               if (z1 .GT. z2 ) then
                  z = z2
                  dz = dz2
               else
                  z = z1
                  dz = dz1
               endif
            else
*              Use inflection starting point if the alternative is to 
*              the right of kmax
               z = z2
               dz = dz2
            endif
         endif
      else
*        The case 1 < p < 2
         z = z2
         dz = dz2
      endif

*     Now solve for kmax
*     Need bounds on either side
      if ( dz .GT. 0.0d00 ) then
*        z  is chosen to be on the left of  kmax
         zlo = z
         zhi = z + 10.d00
         
         flo = dk(p, phi, y, zlo )
         fhi = dk(p, phi, y, zhi )
         allok = 1
      
 565     if ( ( allok .EQ. 1 ) .AND.
     &        ( fhi .GT. 0.d00 ) ) then

            zlo = zhi
            zhi = ( 1.1d00 * zhi ) + 1.0d00

            flo = fhi
            fhi = dk(p, phi, y, zhi )

            if ( zhi .GT. largest/10.d00 ) allok = 0
*              This may not be a problem; it may be that
*              it is so far out no-one really cares!  So
*              let it go, and if problems emerge, they will
*              (hopefully!) be picked up elsewhere.

            goto 565
         endif
      else
         zlo = z / 2.0d00
*        As z is chosen to be on the right of  z

         zhi = z
         flo = dk(p, phi, y, zlo )
         fhi = dk(p, phi, y, zhi )
         allok = 1
      
 566     if ( ( allok .EQ. 1 ) .AND.
     &        ( flo .LT. 0.d00 ) ) then

            zhi = zlo
            zlo = zlo / 2.0d00
            fhi = flo
            flo = dk(p, phi, y, zlo )

            goto 566
         endif
      endif

*     We have the bounds of  zlo  and  zhi
*     To find a useful starting point now, 
*     use linear interpolation:
      z = zlo - flo * (zhi - zlo ) / ( fhi - flo )

*     At last, find the value of tmax (ie. z) and co
      z = sfzero( p, phi, y, zlo, zhi, z,
     &           dk, kddasht, ier )

      tmax = z
      kmax = imgcgf( p, phi, y, tmax )
      
      if ( kmax .LT. 0.0d00 ) then
         kmax = abs(kmax)
         mmax = lrgint
         allok = 0
      else
         dpmmax = ( kmax / pi ) - 0.5d00
         
         if ( dpmmax .GT. dble( lrgint ) ) then
            mmax = lrgint
            allok = 0
         else
            mmax = int( dpmmax )
         endif
      endif

*     If no convergence for tmax, we can probably take it if it is 
*     large enough.  In general, say, the maximum  m  that is used is 5000,
*     when an error is issued.  So if we can make do with m < 5000,
*     take the value as being good enough and continue.

      ier = 0

      return
      end

****************************************************************
*******************************************************************
      double precision function imgdcgf(p, phi, x)
*     Calculates the imag part of the derivative of cgf

      double precision  x, top, bottom, alpha, p, phi, 
     &                  psi, logb

      psi = datan( (1.0d00 - p) * x * phi )
      alpha = 1.0d00 /( 1.0d00 - p )
      top = dcos( psi * alpha )

      logb = alpha  * log( dcos( psi ) )
      bottom = dexp( logb )
*     Appear to need thi fix to get it to work...?      
*        NOTE:  dcos(psi) > 0, so it should work without logs

      imgdcgf = top/bottom

      return
      end

*******************************************************************
*******************************************************************
      double precision function dk(p, phi, y, x)
*     Evaluates the derivative of the  k  function

      double precision  x, p, phi, y

      mu = 1.0d00
      call kdasht( p, phi, y, mu, x, tmp)
      dk = tmp
      return
      end

**************************************************************************
**************************************************************************
      double precision function zerodfn( p, phi, y, x )
*     This evaluates the derivative of the function we are trying to
*     find the zeros of in finding the conditional density.

*     IN:  p, phi, y, mu, x
*     OUT: zerodfn

      double precision  p, mu, phi, y, rl, im, calclambda,
     &                  lambda, x, drl, dim

      mu = 1.0d00
      call kt(p, phi, y, mu, x, rl, im)
      call kdasht(p, phi, y, mu, x, drl, dim)

      lambda = calclambda( p, phi, mu )

      zerodfn = dexp( rl ) * ( -dim * dsin( im ) ) +
     &          dexp( rl )* drl * dcos( im )+
     &          dexp( -lambda ) * y * dsin( x*y )
*      print *,'ZERODFN: rl, im, drl, dim',rl,im,drl,dim

      return
      end

**************************************************************************
****************************************************************
      double precision function g(p, phi, y, x)
*     A function to be numerically integrated, for use with
*     Gauss quadrature method, y > 1.

      double precision  x, p, phi, y, rl, im, imkdash,
     &                  mu, imgdcgf

      mu = 1.0d00
      call kt(p, phi, y, mu, x, rl, im)
      imkdash = imgdcgf(p, phi, x)
      
! *     Compute theta
!       theta = ( mu ** (1.0d00 - p) - 1.0d00 ) / 
!      &         ( 1.0d00 - p ) )
! 
! *     Compute kappa(theta)
!       if ( abs( 2.0d00 - p ) .LT. 1.0d-06 ) 
! *        Use one more term in Taylor series expansion
!          kappa = log( mu ) + 
!      &      (2.0d00 - p) / 2.0d00 * ( ( log(mu) ) ^ 2.0d00 )
!       
!       else
!          mu ** (2.0d00-p) - 1.0d00 ) / ( 2.0d00 - p )
!       endif
! 
!       D = y * theta - kappa

      g = dexp( rl ) / imkdash
      g = dexp( rl ) / ( imkdash - y )
***Changed 03 Dec 1999...look at the definition of imgdcgf.
*  Of couse, if thsi work, CHANGE THE FUNCTION IMGCGF AND ALL CALLS TO IT

      return
      end

*******************************************************************
*******************************************************************
      double precision function intim( p, phi, y, x, m )
*     Computes  Im( k(t) ) - pi/2 - m*pi
*     for finding zeros of the imginary part of he integrand for given  m

      double precision  pi, x, p, phi, y, im, rl, mu
      integer  m

      mu = 1.0d00
      pi = dacos( -1.0d00 )
      call kt(  p, phi, y, mu, x, rl, im )
      intim = -pi/2.0d00 - dble(m) * pi + im

      return
      end

******************************************************************
*******************************************************************
      subroutine findsp( p, mu, phi, y, lowerb, upperb, 
     &                   flo, fhi )

*     Determines a lower and upper bound for the first zero
*     in the density, conditional trick.
***
*     IN:  p, mu, phi, y
*     OUT: lowerb, upperb, flo, fhi
***
      double precision p, mu, phi, y, lowerb, upperb, pi,
     &          calclambda, lambda, t, told, f1, f2,
     &          zerofn, t3, f3, othzero, rl, im,
     &          wt1, wt2, tstep, flo, fhi
      external calclambda

      pi = dacos( -1.0d00 )

***SURELY can improve when y is small, when t gets large.  This
*  routine takes a long time to find the zero.

*     We look at two points of interest.  We have the integrand as
*        dexp( rek ) * dcos( imk ) - dexp( -lambda ) * dcos(t*y)  
*
*     For possible zeros, we look at the first zero of each of the
*     two parts of the function.

*     FIRST, the first zero of the second part of the function, because
*     that's easy:   dexp( -lambda ) * dcos(t*y)

***NECESSARY??? (see later also)
      told = 1.0d-02

      mu = 1.0d00
      t = pi / y
      f2 = zerofn( p, phi, y, t )
*      print *,'First option: ',t, f2

*     SECOND, we now examine the first part of the function for it's
*     first zero.  This is much trickier.

*     We look at the sign of the zeroing function at the turning
*     points and the zeros of the  dsin(t*y)  part.

      f1 = zerofn( p, phi, y, told )
      f2 = zerofn( p, phi, y, t )
*      print *,'f1, zerodfn(0) are: ',f1, zerodfn(p,phi,y,mu,0.0d00)
*AAAAAAAAAAAAAAA
***IT looks to me like  told  is here just to establish the sign of the
*  derivative as we start.  So can we just evaluate  zerodfn() at x=0???

*Well, it appears that  f1  and  zerodfn(..., 0) gives the same sign.
*But I'm not sure it reads any easier.
      call kt(p, phi, y, mu, t, rl, im)
      wt1 = dexp( rl )

      t3 = othzero( p, phi, y )

*     So our two candidate values so far are   t  and  t3.

      f3 = zerofn( p, phi, y, t3 )

      lambda = calclambda(p, phi, mu)
      wt2 = dexp( -lambda )

*     Now we pick the smaller.  It may not be brilliant, but it the only
*     way to guarantee we won't miss the first zero.  The larger one may
*     give a smaller function value, but may be waaaayyy past the first zero.

*      print *,'CANDIDATES:  ', t, t3
      t = min( t, t3 )
      f2 = zerofn( p, phi, y, t )
      tstep = 0.2d00 * t

 100  if ( ( f1 * f2 ) .GT. 0.0d00 ) then

         told = t
         t = told + tstep

         f1 = f2
         f2 = zerofn( p, phi, y, t )
*      print *,'f1*f2 = ',f1*f2

          goto 100

      endif

      lowerb = told
      upperb = t

      flo = f1
      fhi = f2

      return
      end

*******************************************************************
*******************************************************************
      double precision function zerofn( p, phi, y, x )

*     This evaluates the function we are trying to find the zeros
*     of in finding the conditional density.

*     IN:  p, phi, y, x
*     OUT: zerofn
      double precision  p, mu, phi, y, rl, im, calclambda,
     &                  lambda, x

*     When calculating the density, we always have mu=1
      mu = 1.0d00
      call kt(p, phi, y, mu, x, rl, im)
      
      lambda = calclambda( p, phi, mu )

       zerofn = dexp( rl ) * dcos( im )  -
     &         dexp( -lambda ) * dcos( x*y )

*      print *,'p, phi, y, mu, lambda, im, rl, fn1, fn2'
*      print *,p, phi, y, mu, lambda, im, rl, dexp( rl ) * dsin( im ),
*     &        dexp( -lambda ) * dsin( x*y ), zerofn

      return
      end

*****************************************************************


