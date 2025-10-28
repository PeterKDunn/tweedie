
      SUBROUTINE twcomputation(p, phi, y, mu,
     &                 funvalue, exitstatus, relerr, its )

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF.
*
*     IN:   p, phi, y, mu
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION p, phi, y, funvalue, mu
      DOUBLE PRECISION relerr, aimrerr
      DOUBLE PRECISION lambda, Cp, Cy, Cmu, Cphi
      INTEGER exitstatus
      INTEGER its, m
      LOGICAL  pSmall, verbose
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
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
*    funvalue : the value of the function at the given value of  x
*    bound    : The bound using Chebyshev theorm.
*    exitstatus:  1  if relative error is smaller than wished (aimrerr)
*                -1  if not, but the absolute error is less than aimrerr
*               -10  if neither rel or abs error any good
*    Exact zeros acceleration algorithms is used only

*     Defaults
      verbose = .FALSE.
      exitstatus = 1
      relerr = 0.0d0
      its = 0
      
*     Create logical: psmall = TRUE means 1 < p < 2
      psmall = .FALSE.
      IF ( (p .GT. 1.0d00 ) .AND. (p .LT. 2.0d00) ) psmall = .TRUE.
      
*     Compute lambda: Returns lambda = 0.0d0 for p > 2
      CALL findLambda(lambda)

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
      aimrerr = 1.0d-10

      write(*,*) "** Computing for y: ", y
      write(*,*) "**              mu: ", mu
      write(*,*) "**               p: ", p
      write(*,*) "**              phi: ", phi

      IF ( psmall ) THEN
        write(*,*) "About to call DFsmallp from twcomputation"
        CALL DFsmallp(funvalue, exitstatus, relerr, verbose)
      ELSE
        write(*,*) "About to call DFbigp from twcomputation"
        CALL DFbigp(funvalue, exitstatus, relerr, verbose)
      ENDIF
      
*     Fix based on machine accuracy
      write(*,*) "REINSTATE fixes for machine accuracy temp off"
      IF (psmall) THEN
*        IF (funvalue .LT. exp(-lambda) ) funvalue = exp(-lambda)
      ELSE
*        IF (funvalue .LT. 0.0d00) funvalue = 0.0d00
      ENDIF
      write(*,*) "IN twcomputation: funvalue, exitstatus, relerr"
      write(*,*) funvalue, exitstatus, relerr


      RETURN
      END

