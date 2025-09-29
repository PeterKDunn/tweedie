***
*     Author:         Peter Dunn
*     Creation date:  15 February 1999
*     Last revision:  25 Septembeer 2025
*
******************************************************************
** DF-only routines
******************************************************************

      SUBROUTINE DFbigp(funvalue, exitstatus, relerr, its, exact)

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF: p > 2
*
*     IN:   p, phi, y, mu, exacti
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION funvalue, pi
      DOUBLE PRECISION calclambda, resulta, result
      DOUBLE PRECISION result0, relerr, aimrerr
      DOUBLE PRECISION lambda, Cp, Cy, Cmu, Cphi
      DOUBLE PRECISION zeroL, zeroR, zero, kmax, tmax
      DOUBLE PRECISION zeroStartPoint, startPointKmax
      DOUBLE PRECISION omegaInflect
      INTEGER ier, maxit, iteratn, exitstatus
      INTEGER its, exacti, verbose, mfirst, mmax, m, mOld, mNew
      LOGICAL  psmall, exact, stopIterating, convergence
      LOGICAL leftSide
      COMMON /params/ Cp, Cy, Cmu, Cphi, aimrerr
      COMMON /mparam/ m 

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


      write(*,*) " FOR p > 2"
      pi = 4.0d0 * DATAN(1.0d0)

*     FIND kmax, tmax, mmax
      IF (Cy. GT. Cmu) THEN
        write(*,*) "** y >= mu"
        kmax = 0.0d00
        tmax = 0.0d00
        mmax = 0
        mfirst = -1
      write(*,*) "** Found(a): kmax =", kmax
      write(*,*) "             tmax =", tmax
      write(*,*) "             mmax =", mmax
        zeroStartPoint = pi / Cy
      ELSE
        write(*,*) "** y < mu"
        omegaInflect = (pi/2.0d00) * (1.0d00 - Cp)/
     &                 ((2.0d00 * Cp) - 1.0d00)
        startPointKmax = (Cmu ** (1.0d00 - Cp)) / 
     &                   (Cphi * (1.0d00 - Cp)) * 
     &                   DTAN(omegaInflect)
        write(*,*) "Starting pt for finding kmax: ", startPointKmax
        CALL findKmax(kmax, tmax, mmax, mfirst, startPointKmax)
      write(*,*) "** Found(b): kmax =", kmax
      write(*,*) "             tmax =", tmax
      write(*,*) "             mmax =", mmax
      
        IF ( mmax .EQ. 0) THEN
          mfirst = 0
          zeroStartPoint = tmax + pi/Cy
        ELSE
          mfirst = 1
          zeroStartPoint = pi / (Cmu - Cy)
        ENDIF        
      ENDIF
      write(*,*) "           mfirst =", mfirst
      write(*,*) "             StPt =", zeroStartPoint
      write(*,*) "--- (Deal with returned errors, non-convergence)"
      
*     INTEGRATION
      IF (exact) THEN
        write(*,*) "  ----------------------------------"
        write(*,*) "  - Integrating using the EXACT zeroes:"

*       INITAL INTEGRATION
*       For example, even when y >= mu, and Im k(t) heads down, the first
*       regions is not between two zeros. Also true for y < mu.
*
*       IN addition, when y < mu, should we delay integrating until after tmax???
*       Or, accelerate on LEFT of tmax, and RIGHT of tmax, separateky (after 
*       ignoring that first region???

*       INTEGRATE FIRST REGION
        write(*,*) "--- Start integrating FIRST region ---"
        zeroL = (pi / Cy ) * 0.25d00
        zeroR = (pi / Cy ) * 25.0d00
        write(*,*) "DFbigp: zeroL", zeroL
        write(*,*) "DFbigp: zeroR", zeroR
        write(*,*) "DFbigp: zeroSP", zeroStartPoint
*       Now find the right-side zero
        m = mfirst
        CALL findExactZeros(zeroL, zeroR, zeroStartPoint, zero)
        write(*,*) " DFBIGP: FOUND FIRST EXACT ZERO SUCCESSFULLY"        
        zeroL =  0.0d00
        zeroR = zero

*        CALL integrateRegions(zeroL, zeroR)
        write(*,*) "--- DFBigp: Initial integration"
        write(*,*) "    between ", zeroL, "and", zeroR

*       INTEGRATE NON-FIRST REGIONS  
        convergence = .FALSE.
        write(*,*) "--- Standard integration"
 12     IF ( .NOT.(convergence)) THEN
          mOld = m
          m = m - 1
          write(*,*) "  ----------------------------------"
          write(*,*) "    Iterating: m = ",m
          zeroStartPoint = zero
          zeroL = zero
          zeroR = zeroStartPoint * 2.0d00

          CALL findExactZeros(zeroL, zeroR, zeroStartPoint, zero)
          zeroR = zero
          
          write(*,*) "--- DFBigp: integration"
          write(*,*) "    between ", zeroL, "and", zeroR
          
*          CALL integrateRegions()
*          CALLaccelerateExact()
          IF (m .EQ. -4) convergence = .TRUE.
          write(*,*) "--- ABOUT TO RELOOP: convergence", convergence
          GOTO 12
        ENDIF
      ELSE
        write(*,*) "  - DFbifp: for APPROX zeroes:", Cp
        IF ( .NOT.(convergence)) THEN
        write(*,*) "  - Computing for p > 2:", Cp
*           CALL findApproxZeros()
*           CALL integrateRegions()
*           CALL accelerateApprox()
        ENDIF
      ENDIF
      
      RETURN
      END

**************************************************************************

