      SUBROUTINE DFsmallp(funvalue, exitstatus, relerr, exact)

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF: 1 < p < 2
*
*     IN:   p, phi, y, mu, exacti
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION funvalue, pi, zero, zeroL, zeroR
      DOUBLE PRECISION aimrerr, relerr, tmax, kmax, startPoint
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      INTEGER mfirst, m, mOld, mNew, exitstatus, mmax
      LOGICAL  exact, convergence
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

      write(*,*) " FOR 1 < p < 2"
      pi = 4.0d0 * DATAN(1.0d0)
      exitstatus = 0
      relerr = 1.0d00
      convergence = .FALSE.

*     FIND kmax, tmax, mmax
      CALL findKmax(kmax, tmax, mmax, mfirst, startPoint)
      write(*,*) "** Found: mfirst =", mfirst
      write(*,*) "          kmax =", kmax
      write(*,*) "          tmax =", tmax
      write(*,*) "          mmax =", mmax
      
      m = mfirst
      write(*,*) "--- (Deal with returned errors, non-convergence)"

*     INTEGRATION
      if (exact) then
        write(*,*) "  - Computing for EXACT zeroes:"

*       INITAL INTEGRATION
*       For example, even when y >= mu, and Im k(t) heads down, the first
*       regions is not between two zeros. Also true for y < mu.
*
*       IN addition, when y < mu, should we delay integrating until after tmax???
*       Or, accelerate on LEFT of tmax, and RIGHT of tmax, separateky (after 
*       ignoring that first region???


*       INTEGRATE FIRST REGION
        write(*,*) "--- Start integrating ---"
        zeroL =  (pi / Cy ) / 2.0d00
        zeroR = (pi / Cy ) * 2.0d00

*       Now find the right-side zero
        CALL findExactZeros(zeroL, zeroR, zero)
        zeroL =  0.0d00
        zeroR = zero
        
*        CALL integrateRegions(zeroL, zeroR)
        write(*,*) "--- Initial integration"
        write(*,*) "    between ", zeroL, "and", zeroR


*       INTEGRATE NON-FIRST REGIONS      
        write(*,*) "--- Standard integration"
        if ( .NOT.(convergence)) then
          mOld = m
          CALL advanceM(mmax, kmax, leftSide, mOld, mNew)
          m = mNew
          
          write(*,*) "    Solving, using m:", m
          
        CALL findExactZeros(zeroL, zeroR, zero)
*          CALL integrateRegions()
*          CALLaccelerateExact()
        endif
      else
        write(*,*) "  - Computing for APPROX zeroes:", Cp
        if ( .NOT.(convergence)) then
        write(*,*) "  - Computing for p > 2:", Cp
*           CALL findApproxZeros()
*           CALL integrateRegions()
*           CALL accelerateApprox()
        endif
      endif
      
      RETURN
      END

**************************************************************************

