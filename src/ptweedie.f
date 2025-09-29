      SUBROUTINE twcdf(p, phi, y, mu, exacti,
     &                 funvalue, exitstatus, relerr, its )

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF.
*
*     IN:   p, phi, y, mu, exacti
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION p, phi, y, funvalue, mu, pi
      DOUBLE PRECISION calclambda, resulta, result
      DOUBLE PRECISION result0, relerr, aimrerr
      DOUBLE PRECISION lambda, Cp, Cy, Cmu, Cphi
      DOUBLE PRECISION zeroL, zeroR, zero
      DOUBLE PRECISION kmax, tmax, startPoint
      INTEGER ier, maxit, iteratn, exitstatus
      INTEGER its, exacti, verbose, mfirst, mmax, m, mOld, mNew
      LOGICAL  psmall, exact, stopIterating, convergence
      LOGICAL leftSide
      COMMON /params/ Cp, Cy, Cmu, Cphi, aimrerr
      COMMON /mparam/ m 

*     Set COMMON block vars
      Cp = p
      Cmu = mu
      Cphi = phi
      Cy = y

* VARIABLES:
*    y        : the value at which the function is to be evaluated
*    x        : the range over which the integral is to be integrated; an
*               internal variable; NOT the value at which the function is to be evaluated
*    lambda   : for 1 < p < 2, P(Y = 0) = exp( -lambda )
*    p        : the index (i.e., variance function is V(mu) = mu ^ p)
*    phi      : the dispersion parameter
*    funvalue : the value of the function at the given value of  x
*    bound    : The bound using Chebyshev theorm.
*    exitstatus:  1  if relative error is smaller than wished (aimrerr)
*                -1  if not, but the absolute error is less than aimrerr
*               -10  if neither rel or abs error any good
*    exact    : 1 if the exact zeros acceleration algorithms is used;
*               0 if the approx zeros acceleration algorithm is used.

      pi = 4.0d0 * DATAN(1.0d0)

*     Defaults
      verbose = 0
      exitstatus = 1
      relerr = 0.0d00
      its = 0
      
*     Create logical: exact = .TRUE. means use exact zeros in acceleration
      exact = .TRUE.
      IF (exacti. EQ. 0) exact = .FALSE.

*     Create logical: psmall = TRUE means 1 < p < 2
      lambda = 0.0d00
      psmall = .FALSE.
      IF ( (p .GT. 1.0d00 ) .AND. (p .LT. 2.0d00) ) THEN
        psmall = .TRUE.
        CALL findLambda(lambda, p, mu, phi)
      ENDIF

*     SPECIAL CASE: if y < 0, return 0
      IF ( y .LT. 0.0d00 ) then
        funvalue = 0.0d00
        RETURN
      ENDIF

*     SPECIAL CASE: if 1 < p < 2, Pr(Y = 0) = exp( -lambda )
      IF ( psmall .AND. (y .EQ. 0.0d00 ) ) THEN
        funvalue = DEXP( -lambda )
        RETURN
      ENDIF

*     SET ACCURACY REQUIREMENTS
      maxit = 25
      write(*,*) ">>>> Consider changing maxit here"
      aimrerr = 1.0d-10

*     set other parameters
      iteratn = 0
      ier = 0
      resulta = 0.0d00
      result  = 0.0d00
      result0 = 0.5d00
      convergence = .FALSE.
      stopIterating = .TRUE.

      write(*,*) "** Computing for y ", y
      write(*,*) "**               mu ", mu
      write(*,*) "**               p ", p
      write(*,*) "**               phi ", phi

      IF ( psmall ) THEN
        CALL DFsmallp(funvalue, exitstatus, relerr, its, exact)
      ELSE
        CALL DFbigp(funvalue, exitstatus, relerr, its, exact)
      ENDIF


    
      RETURN
      END

