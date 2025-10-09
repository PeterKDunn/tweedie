
      SUBROUTINE DFbigp(funvalue, exitstatus, relerr, exacti, verbose)

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF: p > 2
*
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION funvalue, pi, sum
      DOUBLE PRECISION relerr, aimrerr, epsilon
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi, areaT, Wold2
      DOUBLE PRECISION zeroL, zeroR, zero, kmax, tmax
      DOUBLE PRECISION zeroStartPoint, startTKmax, wvec(300)
      DOUBLE PRECISION West, Wold, area0, xvec(300)
      DOUBLE PRECISION area1, areaA, psi
      DOUBLE PRECISION zeroBoundL, zeroBoundR
      DOUBLE PRECISION DFintegrand, findKmaxSP
      DOUBLE PRECISION Mmatrix(2, 200), Nmatrix(2, 200)
      EXTERNAL DFintegrand, findKmaxSP
      INTEGER exitstatus, itsAcceleration, itsPreAcc
      INTEGER mfirst, mmax, m, mOld, accMax, exacti
      LOGICAL exact, convergence, leftOfMax, verbose
      LOGICAL stopPreAccelerate, pSmall, flip
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      COMMON /mparam/ m 

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
*    exact    : 1 if the exact zeros acceleration algorithms is used;
*               0 if the approx zeros acceleration algorithm is used.


      write(*,*) " FOR p > 2"
      
      pi = 4.0d0 * DATAN(1.0d0)
      exitstatus = 0
      relerr = 1.0d00
      epsilon = 1.0d-16
      aimrerr = 1.0d-14
      convergence = .FALSE.
      exact = .TRUE.
      exacti = 1
      verbose = .TRUE.
      
*     FIND kmax, tmax, mmax
      IF (Cy. GE. Cmu) THEN
*     ************** y > MU   **************
        write(*,*) "** y >= mu"
        kmax = 0.0d00
        tmax = 0.0d00
        mmax = 0
        mfirst = -1
        mOld = 0
        write(*,*) "** Im k(t) heads down immediately"
      
        zeroStartPoint = pi / Cy
        leftOfMax = .FALSE.
      ELSE
*     ************** y < MU   **************
        write(*,*) "** y < mu"
        
        startTKMax = findKmaxSP()

        write(*,*) "Starting t for finding kmax: ", startTKmax
        CALL findKmax(kmax, tmax, mmax, mfirst, startTKmax)
        
        write(*,*) "** Found(b): kmax =", kmax
        write(*,*) "             tmax =", tmax
        write(*,*) "             mmax =", mmax

        leftOfMax = .TRUE.
        IF ( mmax .EQ. 0) THEN
          mfirst = 0
          mOld = 0
          zeroStartPoint = tmax + pi/Cy
          leftOfMax = .FALSE.
        ELSE
          mfirst = 1
          mOld = 0
          zeroStartPoint = pi / (Cmu - Cy)
          mOld = m

          CALL advanceM(mmax, m, mOld, leftOfMax, flip)

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
*      write(*,*) "    --- Find right-side zero for m:", mfirst
*      write(*,*)" ALREADY HAVE: ", zeroStartPoint
      zeroBoundL = 0.0d00
      zeroBoundR = zeroStartPoint * 2.0d00

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

      itsPreAcc = 0
      IF (mfirst .EQ. -1 ) THEN
*       Accelerate immediately; 'no pre-acceleration' area
        itsPreAcc = itsPreAcc + 1
        write(*,*) "  > Not using pre-acceleration area"
        area1 = 0.0d00

        mOld = m

        CALL advanceM(mmax, m, mOld, leftOfMax, flip)

      ELSE
*       Find some areas BEFORE accelerating
        area1 = 0.0d00

        mOld = m

        CALL advanceM(mmax, m, mOld, leftOfMax, flip)

        stopPreAccelerate = .FALSE.
 115    IF ( .NOT.(stopPreAccelerate) ) THEN
          itsPreAcc = itsPreAcc + 1

          IF (leftOfMax ) THEN
             zeroBoundL = zeroR
             zeroBoundR = zeroR * 10.00d00
          ELSE
            zeroBoundL = tmax 
            zeroBoundR = zeroR * 20.0d00
          ENDIF
          zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0d00
*       write(*,*) "--> Start pt", zeroStartPoint
*       write(*,*) "--> BoundsLt", zeroBoundL
*       write(*,*) "--> BoundsR", zeroBoundR
          zeroL = zeroR
          CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                        zeroStartPoint, zero)

          zeroR = zero
          write(*,*) "--- Integrate (m = ", m, ") between " 
  
          CALL gaussq( DFintegrand, sum, zeroL, zeroR)
          area1 = area1 + sum
          write(*,*) zeroL, "and ", zeroR, ": ", sum

*         STOP condition for pre-acceleration.
*         Not sure about this...
*          if ( m .EQ. (mmax - 1) ) stopPreAccelerate = .TRUE.
          if ( itsPreAcc .GE. 2) stopPreAccelerate = .TRUE.
          CALL advanceM(mmax, m, mOld, leftOfMax, flip)

          GOTO 115
        ENDIF
      ENDIF

*      write(*,*) "Finished pre-acc; areas"
*      write(*,*) "SUMMARY (before accelerating):"
*      write(*,*) "  Area0 ", area0
*      write(*,*) "  Area1 ", area1

*     3. INTEGRATE: the ACCELERATION regions: areaA
      write(*,*) "*******************************" 
      write(*,*) "3. INTEGRATE: the ACCELERATION"
      
      Wold = 0.0d00
      Wold2 = 1.0d00

      IF (exact) THEN

        itsAcceleration = 0
        areaA = 0.0d00
        convergence = .FALSE.

        xvec(1) = zeroR
*       This will be the very first, left-most value of t used, the left-most
*       value of  t  used in the acceleration (the previous regions *right* value) 

 12     IF ( .NOT.(convergence)) THEN
          write(*,*) "  --- Next tail region"

          itsAcceleration = itsAcceleration + 1
*         itsAcceleration = 1 means this is the first area found
*         under the acceleration regime

          IF (leftOfMax ) THEN
            zeroStartPoint = zeroR
            zeroL = zeroR 
            zeroR = zeroR * 20.0d00
          ELSE
            IF (flip) THEN
*             FLIPPING to other side of tmax

              zeroStartPoint = tmax + ( tmax - zero)
*             That is, start of the other side of tmax            
              zeroL = zero
              zeroR = zeroStartPoint * 20.0d00
            ELSE
              zeroStartPoint = zeroR
              zeroL = zeroR
              zeroR = zeroR * 10.0d00
            ENDIF
          ENDIF
*  write(*,*) "that factor of 20: depeds on slope!"
* write(*,*) "Flatter? Larger multiplier"
*        write(*,*) "Steeper? Smaller multiplier"

          CALL findExactZeros(zeroL, zeroR, 
     &                        zeroStartPoint, zero)
          IF (leftOfMax) THEN
            zeroR = zero
          ELSE
            zeroR = zero
          ENDIF

          xvec(itsAcceleration + 1) = zeroR

          write(*,*) "  - Integrate (m = ", m, "):", zeroL, zeroR

          CALL gaussq( DFintegrand, psi, zeroL, zeroR)
*         psi: area of the latest region
          wvec(itsAcceleration) = psi
          write(*,*) "  - Area between zeros is:", psi

          accMax = 40
          Wold2 = Wold
          Wold = West
          CALL accelerateNEW(xvec, wvec, itsAcceleration, 
     &     Mmatrix, Nmatrix, West)
*          W is the best guess of the convergent integration
          write(*,*) "iteration", itsAcceleration, ":", West
          write(*,*) "  - Estimate of tail area:", West
          write(*,*) "--------------------------------"

*         Check for convergence
          relerr = (DABS( West - Wold ) + DABS( West - Wold2 ) ) /
     &                     (DABS(West) + epsilon )
          IF (relerr .LT. aimrerr ) THEN 
            write(*,*) "  Relerr is", relerr
            convergence = .TRUE.
          ENDIF
*        IF (m .EQ. -10) convergence = .TRUE.

          mOld = m
          CALL advanceM(mmax, m, mOld, leftOfMax, flip)

          GOTO 12
        ENDIF
      ELSE
        IF ( .NOT.(convergence)) THEN
        write(*,*) "  - Computing for p > 2:", Cp
*           CALL findApproxZeros()
*           CALL integrateRegions()
*           CALL accelerate()
        ENDIF
      ENDIF
      
      areaA = West
      areaT = area0 + area1 + areaA
      write(*,*) "SUMMARY:"
      write(*,*) "  Area0 ", area0
      write(*,*) "  Area1 ", area1, "(", itsPreAcc, "regions)"
      write(*,*) "  AreaA ", areaA, "(", itsAcceleration, " its)"
      write(*,*) "  TOTAL ", areaT
      
      
*** WHAT DTO DO with relerrr? Might have three rel eerrors: from initila, pre-acc, acc?
*** Take largest of the three? ADD?
      write(*,*) "FIX rel err: |A|.relA + ... + |C|.relC/|A+B+C|"
      
*     We have the value of the integral in the CDF calculation. 
*     So now work out the CDF
      funvalue = (-1.0d00/pi) * areaT + 0.5d00
      write(*,*) "FINAL AREA: The cdf value is", funvalue
      write(*,*) "DFbigp: funvalue, exitstatus, relerr, exacti"
      write(*,*) funvalue, exitstatus, relerr, exacti

      RETURN
      END

**************************************************************************

