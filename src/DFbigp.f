      SUBROUTINE DFbigp(funvalue, exitstatus, relerr, exact)

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF: p > 2
*
*     IN:   p, phi, y, mu, exacti
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION funvalue, pi, sum
      DOUBLE PRECISION relerr, aimrerr, epsilon
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi, areaT, Wold2
      DOUBLE PRECISION zeroL, zeroR, zero, kmax, tmax
      DOUBLE PRECISION zeroStartPoint, startTKmax, wvec(200)
      DOUBLE PRECISION omegaInflect, DFintegrand, West, Wold
      DOUBLE PRECISION area0, area1, areaA, psi, xvec(200)
      DOUBLE PRECISION zeroBoundL, zeroBoundR
      EXTERNAL DFintegrand
      INTEGER exitstatus, itsAcceleration
      INTEGER mfirst, mmax, m, mOld, accMax
      LOGICAL exact, convergence, leftOfMax
      LOGICAL stopPreAccelerate
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
      exitstatus = 0
      relerr = 1.0d00
      epsilon = 1.0d-16
      aimrerr = 1.0d-14
      convergence = .FALSE.
      
*     FIND kmax, tmax, mmax
      IF (Cy. GE. Cmu) THEN
*     ************** y > MU   **************
        write(*,*) "** y >= mu"
        kmax = 0.0d00
        tmax = 0.0d00
        mmax = 0
        mfirst = -1
      write(*,*) "** Found(a): kmax =", kmax
      write(*,*) "             tmax =", tmax
      write(*,*) "             mmax =", mmax
        zeroStartPoint = pi / Cy
        leftOfMax = .FALSE.
      ELSE
*     ************** y < MU   **************
        write(*,*) "** y < mu"
        omegaInflect = (pi/2.0d00) * (1.0d00 - Cp)/
     &                 ((2.0d00 * Cp) - 1.0d00)
        startTKmax = (Cmu ** (1.0d00 - Cp)) / 
     &                   (Cphi * (1.0d00 - Cp)) * 
     &                   DTAN(omegaInflect)

        write(*,*) "Starting t for finding kmax: ", startTKmax
        CALL findKmax(kmax, tmax, mmax, mfirst, startTKmax)
        write(*,*) "** Found(b): kmax =", kmax
        write(*,*) "             tmax =", tmax
        write(*,*) "             mmax =", mmax

        leftOfMax = .TRUE.
        IF ( mmax .EQ. 0) THEN
          mfirst = 0
          zeroStartPoint = tmax + pi/Cy
          leftOfMax = .FALSE.
        ELSE
          mfirst = 1
          zeroStartPoint = pi / (Cmu - Cy)
          mOld = m
          CALL advanceM(mmax, m, mOld, leftOfMax)
        ENDIF        
      ENDIF
  
      write(*,*) "           mfirst =", mfirst
      write(*,*) "             StPt =", zeroStartPoint
      write(*,*) "--- (Deal with returned errors, non-convergence)"



*     INTEGRATION
*     There are three integration regions:
*
*     1. The *initial* area, which is not between zeros of Im{k(t)}: area0
*
*     2. The initial area *before* Sidi acceleration is invoked: area1
*        (For instance, wait until Im{k(t)} is on the downturn.)
*
*     3. The area thereafter, upon which Sidi acceleration is
*        applied; the area returned by acceleration is areaA

      area0 = 0.0d00
      area1 = 0.0d00
      areaA = 0.0d00

*     1. INTEGRATE FIRST REGION: area0
      write(*,*) "*******************************" 
      write(*,*) "1. INTEGRATE: the INITIAL region"
      write(*,*) "    --- Find the right-side zero"
      write(*,*) "Using m = ", mfirst
      zeroBoundL = 0.0d00
      zeroBoundR = (pi / Cy ) * 2.0d00

*     Now find the right-side zero
      m = mfirst
      CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                    zeroStartPoint, zero)
      zeroL =  0.0d00
      zeroR = zero

      write(*,*) "  - Between ", zeroL, zeroR
      CALL gaussq( DFintegrand, area0, zeroL, zeroR)
      write(*,*) "  - Initial area is", area0
      

      write(*,*) "*******************************" 
*     2. INTEGRATE: the PRE-ACCELERATION regions: area1
      write(*,*) "2. INTEGRATE: the PRE-ACCELERATION regions"
*     When p > 2, things seem well-behaved most of the time, so 
*     we declare  area1  to be up to m = mmax - 1 (i.e., just
*     after the downturn)

      IF (mfirst .EQ. -1 ) THEN
*       Accelerate immediately; 'no pre-acceleration' area
        write(*,*) "  > Not using pre-acceleration area"
        area1 = 0.0d00

        mOld = m
        CALL advanceM(mmax, m, mOld, leftOfMax)
      ELSE
*       Find some areas BEFORE accelerating
        area1 = 0.0d00

        mOld = m
        CALL advanceM(mmax, m, mOld, leftOfMax)

        stopPreAccelerate = .FALSE.
 115    IF ( .NOT.(stopPreAccelerate) ) THEN
          IF (leftOfMax ) THEN
             zeroBoundL =  zeroR
             zeroBoundR = tmax
          ELSE
            zeroBoundL = tmax 
            zeroBoundR = zeroR * 20.0d00
          ENDIF
          zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0d00

          zeroL = zeroR
          CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                        zeroStartPoint, zero)

          zeroR = zero
          write(*,*) "--- Integrate (m = ", m, ") between " 
          write(*,*) zeroL, "and ", zeroR
  
          CALL gaussq( DFintegrand, sum, zeroL, zeroR)
          area1 = area1 + sum
          write(*,*) "... giving ", sum
          
          if ( m .EQ. -1 ) stopPreAccelerate = .TRUE.
          CALL advanceM(mmax, m, mOld, leftOfMax)

          GOTO 115
        ENDIF
      ENDIF

      write(*,*) "Finished pre-acc; areas"
      write(*,*) "SUMMARY (before accelerating):"
      write(*,*) "  Area0 ", area0
      write(*,*) "  Area1 ", area1


*     3. INTEGRATE: the ACCELERATION regions: areaA
      write(*,*) "*******************************" 
      write(*,*) "3. INTEGRATE: the ACCELERATION"
      
      Wold = 0.0d00
      Wold2 = 1.0d00

      IF (exact) THEN
        write(*,*) "----------------------------------"
        write(*,*) "  - Integrating using the EXACT zeroes:"

        itsAcceleration = 0
        areaA = 0.0d00
        convergence = .FALSE.

        xvec(1) = zeroR
*       This will be the very first, left-most value of t used, the left-most
*       value of  t  used in the acceleration (the previous regions *right* value) 

 12     IF ( .NOT.(convergence)) THEN
          write(*,*) "            "
          write(*,*) "INTEGRATE: all subsequent regions"

          itsAcceleration = itsAcceleration + 1
*         itsAcceleration = 1 means this is the first area found
*         under the acceleration regime


          zeroStartPoint = zeroR
          zeroL = zeroR 
          zeroR = zeroR * 20.0d00

          CALL findExactZeros(zeroL, zeroR, 
     &                        zeroStartPoint, zero)
          zeroR = zero
          xvec(itsAcceleration + 1) = zeroR

          write(*,*) "--- Integrate (m = ", m, ") between " 
          write(*,*) zeroL, "and ", zeroR
          
          CALL gaussq( DFintegrand, psi, zeroL, zeroR)
*         psi: area of the latest region
          wvec(itsAcceleration) = psi
          write(*,*) "  - Area between zeros is:", psi

          accMax = 40
          Wold2 = Wold
          Wold = West
          CALL accelerate(xvec, wvec, itsAcceleration, 
     &                    accMax, West)
*          W is the best guess of the convergent integration
          write(*,*) "iteration", itsAcceleration, ":", West
          write(*,*) "  - Estimate of tail area:", West
          write(*,*) "--------------------------------"

*         Check for convergence
         relerr = (DABS( West - Wold ) + DABS( West - Wold2 ) ) /
     &           (DABS(West) + epsilon )
          IF (relerr .LT. aimrerr ) THEN 
            write(*,*) "  Relerr is", relerr
            convergence = .TRUE.
          ENDIF

          mOld = m
          CALL advanceM(mmax, m, mOld, leftOfMax)
* THIS SHOULD JUST BE m = m = 1 BY NOW...



          GOTO 12
        ENDIF
      ELSE
        write(*,*) "  - DFbifp: for APPROX zeroes:", Cp
        IF ( .NOT.(convergence)) THEN
        write(*,*) "  - Computing for p > 2:", Cp
*           CALL findApproxZeros()
*           CALL integrateRegions()
*           CALL accelerate()
        ENDIF
      ENDIF
      
      areaA = West
      areaT = area 0 + area1 + areaA
      write(*,*) "SUMMARY:"
      write(*,*) "  Area0 ", area0
      write(*,*) "  Area1 ", area1
      write(*,*) "  AreaA ", areaA
      write(*,*) "  TOTAL ", areaT
      
*     We have the value of the integral in the CDF calculation. 
*     So now work out the CDF
      funvalue = (-1.0d00/pi) * areaT + 0.5
      write(*,*) "FINAL AREA: The cdf value is", funvalue

      RETURN
      END

**************************************************************************

