***
*     Author:         Peter Dunn
*     Creation date:  15 February 1999
*     Last revision:  25 Septembeer 2025
*
* IN R, print statements for debugging as follows:
*     call dblepr("The value of y is ",-1, y, 1)
*     call intpr("The value of exact is ",-1, exact, 1)
******************************************************************
** Routines common to density and DF
******************************************************************


*******************************************************************
      double precision function calclambda(p, phi, mu )
*     Calculates lambda, such that  P(Y = 0) = exp( -lambda )

      double precision  p,  phi, mu, twomp

      twomp = 2.0d00 - p      
      calclambda = (mu ** twomp) / (phi * twomp)

      return
      end

*******************************************************************
****************************************************************
      subroutine kt(p, phi, y, mu, x, rl, im)
*     Calculates real and imaginary parts of  k(t)  for the DF: k(t) = K(it) - ity
*
*     IN:  p, phi, y, mu, x
*     OUT: rl, im
*
      double precision  p, phi, y, mu, x, rl, im,
     &                  omega, alpha, front, denom

* MAJOR VARIABLES:
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf

      omega = datan( ( 1.0d00 - p ) * x * phi /
     &               ( mu ** (1.0d00 - p) ) )

      alpha = ( 2.0d00 - p ) / ( 1.0d00 - p )
      front = ( mu ** ( 2.0d00 - p ) ) / 
     &         ( phi * ( 2.0d00 - p ) )
      denom = dcos( omega ) ** alpha

      rl = front * ( dcos ( omega * alpha ) / denom - 1)
      im = ( front * dsin ( omega * alpha ) / denom ) - 
     &         x * y

      return
      end

****************************************************************
****************************************************************
      subroutine kdasht(p, phi, y, mu, x, rl, im)
*     Calculates real and imaginary parts of  k'(t)  for the DF: k(t) = K(it) - ity
*
*     IN:  p, phi, y, mu, x
*     OUT: rl, im

      double precision  p, phi, y, mu, x, rl, im,
     &                  omega, alpha, denom

* MAJOR VARIABLES:
*   rl         : the real part of the derivative of the cgf
*   im         : the imaginary part of the derivative of the cgf

      omega = datan( ( 1.0d00 - p ) * x * phi /
     &               ( mu ** (1.0d00 - p) ) )

      alpha = 1.0d00 / ( 1.0d00 - p )
      denom = dcos( omega ) ** alpha

      rl = -mu * ( dsin( omega * alpha ) / denom )
      im =  mu *   dcos( omega * alpha ) / denom - y

      return
      end

****************************************************************
*****************************************************************
      subroutine kddasht( p, phi, mu, x, ddrl, ddim)
*     Calculates real and imaginary parts of  k''(t)  for the DF: k(t) = K(it) - ity

*     IN:   p, phi, mu, x
*     OUT:  ddrl, ddim

      double precision  p, phi, mu, x, ddrl, ddim, 
     &                  denom, alpha, omega

* MAJOR VARIABLES:
*     ddrl      : the real part of the second derivative
*     ddim      : the imaginary part of the second derivative

      omega = datan( ( 1.0d00 - p ) * x * phi /
     &               ( mu ** (1.0d00 - p) ) )

      alpha = p / ( 1.0d00 - p )

      denom = dcos( omega ) ** alpha

      ddrl = -phi * mu**p *    dsin( omega * alpha ) 
     &            / denom
      ddim = -phi * mu**p  * ( dcos( omega * alpha ) 
     &            / denom )

      return
      end

******************************************************************
*******************************************************************
      subroutine accelerate( FF, psi, xvec, mmatrix, nmatrix, w,
     &                    znum, relerr, wold, sumarea, 
     &                    flag, verbose )
*     Accelerates the integration using Sidi's (1988) method. Faster 
*     than Shanks' method for the series with which we are concerned.

      double precision  mmatrix(2, 200), nmatrix(2, 200),
     &          xvec(200), FF, psi, w, relerr, wold(3),
     &          denom, sumarea, abserr, largest, wsave
      integer  i, znum, ell, flag, verbose

*     VARIABLES
*     FF      : The cdf up to x_l
*     psi     : the integral between x_l and x_{l+1}
*     xvec    : contains the x-values used up to  x_{l+1}
*     znum    : the iteration number for the table
*     sumarea : the sum of the areas; rarely, psi = 0.0, which
*               causes things to die.  In this case, use sumarea.
*     flag    : if FLAG=1,  then the limits of the machine are being
*               reached (largest integer).
*     wsave   : the old value of w in case something goes awry
*     w       : the value of the series summation, since the first
*               value of x used.

      ell = znum - 1
      flag = 0
      largest = 1.0d30

      if ( abs(psi) .LT. 1.0d-31 ) then
         w = FF
         relerr = 0.0d00

         return
      else

         mmatrix(2, 1) = FF / psi
         nmatrix(2, 1) = 1.0d00 / psi
         sumarea = sumarea + psi
*         call intpr("In accelerate, verbose = ", -1, verbose, 1)
         
         if ( verbose .EQ. 1 ) then
            call dblepr("    w(x) = ", -1, psi, 1)
            call dblepr("    F(x) = ", -1, FF, 1)
            call dblepr("    M-matrix (2,1) = ", -1, 
     &                  mmatrix(2,1), 1)
            call dblepr("    N-matrix (2,1) = ", -1, 
     &                  nmatrix(2,1), 1)
         endif

*        Add the new information
         flag = 0
         do i = 2, znum
   
            if ( verbose .EQ. 1 ) then
               call intpr("    Adding new info at element ", 
     &                     -1, i, 1)
            endif
            
            denom = 1.0d00 / xvec(znum + 1 - i) - 
     &                        1.0d00 / xvec(znum)
   
            mmatrix(2, i) = ( mmatrix(1, i - 1) - 
     &                        mmatrix(2, i - 1) ) 
     &                      / denom
            nmatrix(2, i) = ( nmatrix(1, i - 1) - 
     &                        nmatrix(2, i - 1) ) 
     &                      / denom

            
            if ( verbose .EQ. 1 ) then
               call dblepr("    demoninator = ", -1, denom, 1)
               call dblepr("    New M-matrix entry = ", -1, 
     &                     mmatrix(2,i), 1)
               call dblepr("    New N-matrix entry = ", -1, 
     &                     nmatrix(2,i), 1)
            endif
            
            if ( (abs(mmatrix(2, i)) .GT. largest) .OR.
     &           (abs(nmatrix(2, i)) .GT. largest) ) then
               flag = 1
            endif
         enddo

         if ( (abs(mmatrix(2,znum)) .GT. largest) .OR.
     &        (abs(nmatrix(2,znum)) .GT. largest) ) then
            flag = 1
            wsave = w
         else
            if ( znum .GT. 1 ) then
               w = mmatrix(2, znum) / nmatrix(2, znum)
               if ( verbose .EQ. 1 ) then
                  call dblepr("    New W value = ", -1, w, 1)
               endif
            endif

            wold(1) = wold(2)
            wold(2) = wold(3)
            wold(3) = w

         endif

*        This is the error used in DQEXT:
*         relerr = abs( w - wold(1) ) + abs( ( w - wold(2) ) )
*        This causes problems here though, and failure to work?????
*        (eg, try it with p=1.5, mu=phi=1, y=0(100)10)?????????????

         if ( ell .GT. 3 ) then
            relerr = abs( w - wold(1) ) + 
     &               abs( ( w - wold(2) ) ) / w
            abserr = abs( wold(3) - wold(2) )
            if ( verbose .EQ. 1 ) then
               call dblepr("    Rel. error estimate = ", 
     &                     -1, relerr, 1)
            endif
         else
            relerr = 1.0d00
         endif

         if ( w .NE. 0.0d00 ) then
            relerr = relerr
         endif

*        Now drop first row
         do  i = 1, znum
            mmatrix(1, i) = mmatrix(2, i)
            nmatrix(1, i) = nmatrix(2, i)
         enddo
      endif

      return
      end
  
