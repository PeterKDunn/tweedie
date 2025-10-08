
      SUBROUTINE DFsmallp(funvalue, exitstatus, relerr, exact)

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF: 1 < p < 2
*
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION funvalue, pi, zero, zeroL, zeroR, sum
      DOUBLE PRECISION aimrerr, relerr, tmax, kmax, f, df
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi, findKmaxSP, startTKMax
      DOUBLE PRECISION zeroStartPoint, area0, area1, areaA
      DOUBLE PRECISION zeroBoundL, zeroBoundR, DFintegrand, psi
      DOUBLE PRECISION Wold, Wold2, areaT, epsilon
      DOUBLE PRECISION West, xvec(300), wvec(300), lambda
      DOUBLE PRECISION Mmatrix(2, 200), Nmatrix(2, 200)
      DOUBLE PRECISION finalTP
      INTEGER mfirst, m, mOld, exitstatus, mmax, n, i
      INTEGER itsPreAcc, accMax, exacti, itsAcceleration
      LOGICAL  exact, convergence, flip, leftOfMax
      LOGICAL pSmall, stopPreAccelerate
      EXTERNAL findKmaxSP, DFintegrand
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

      write(*,*) " FOR 1 < p < 2"
      
*     Initialise the M and N matrices, x and w vectors
      DO i = 1, 200
        Mmatrix(1, i) = 0.0d0
        Mmatrix(2, i) = 0.0d0
        Nmatrix(1, i) = 0.0d0
        Nmatrix(2, i) = 0.0d0
        xvec(i) = 0.d000
        wvec(i) = 0.0d00
      ENDDO
      
*     Set up
      pi = 4.0d0 * DATAN(1.0d0)
      exitstatus = 0
      relerr = 1.0d00
      epsilon = 1.0d-16
      aimrerr = 1.0d-14
      convergence = .FALSE.
      exact = .TRUE.
      exacti = 1


*     FIND kmax, tmax, mmax
      IF (Cy. GE. Cmu) THEN
*     ************** y > MU   **************
*     Im k(t) heads down immediately

        write(*,*) "** y >= mu"
        kmax = 0.0d00
        tmax = 0.0d00
        mmax = 0
        mfirst = -1
        mOld = 0
      
        zeroStartPoint = pi / Cy
        leftOfMax = .FALSE.
      ELSE
*     ************** y < MU   **************
*       HARDER!         
        write(*,*) "** y < mu"
        
        startTKMax = findKmaxSP()

        write(*,*) "Starting t for finding kmax: ", startTKmax
        CALL findKmax(kmax, tmax, mmax, mfirst, startTKmax)
        
        write(*,*) "** Found(b): kmax =", kmax
        write(*,*) "             tmax =", tmax
        write(*,*) "             mmax =", mmax

*        leftOfMax = .TRUE.
*        IF ( mmax .EQ. 0) THEN
*          mfirst = 0
*          mOld = 0
*          zeroStartPoint = tmax + pi/Cy
*          leftOfMax = .FALSE.
*        ELSE
*          mfirst = 1
*          mOld = 0
*          zeroStartPoint = pi / (Cmu - Cy)
*          mOld = m
*
*          CALL advanceM(mmax, m, mOld, leftOfMax, flip)
*        ENDIF
      ENDIF
*      
*      write(*,*) "           mfirst =", mfirst
*      write(*,*) "             StPt =", zeroStartPoint
*      write(*,*) "--- (Deal with returned errors, non-convergence)"


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
      m = mfirst 
      
      write(*,*) "ABOUT TO FINAL TURNING"
      
*     Find the final turning point of Im/Re k, and start accelerating thereafter      
      write(*,*) "IS this about TPs correct??"
      finalTP = 0.0d00
      IF ( Cy .LT. Cmu ) THEN
        CALL findAccelStart(finalTP)
      ENDIF
      
      write(*,*) "START ACCELERATION AFTER: t:", finalTP


      zeroStartPoint = pi / Cy
*      write(*,*) "Start pt for first zero:", zeroStartPoint

      
*     1. INTEGRATE FIRST REGION: area0
      write(*,*) "*******************************" 
      write(*,*) "1. INTEGRATE: the INITIAL region"
*      write(*,*) "    --- Find right-side zero"
*      write(*,*)" ALREADY HAVE: ", zeroStartPoint
      zeroBoundL = zeroStartPoint - 0.25d0 * pi / Cy
      zeroBoundR = zeroStartPoint + 0.25d0 * pi / Cy
*      write(*,*) "    & is between ", zeroBoundL, zeroBoundR

*     Now find the right-side zero
      CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                    zeroStartPoint, zero)
      write(*,*) "ZERO: ", zero
      zeroL =  0.0d00
      zeroR = zero

      CALL gaussq( DFintegrand, area0, zeroL, zeroR)
      write(*,*) "  - Initial area is", area0
      write(*,*) "    between ", zeroL, " and ", zeroR
      

      write(*,*) "*******************************" 
*     2. INTEGRATE: the PRE-ACCELERATION regions: area1
      write(*,*) "2. INTEGRATE: the PRE-ACCELERATION regions"
*     When p > 2, things seem well-behaved most of the time, so 
*     we declare  area1  to be up to m = mmax - 1 (i.e., just
*     after the downturn)

      itsPreAcc = 0
      area1 = 0.0d00
      CALL advanceM(mmax, m, mOld, leftOfMax, flip)

*      IF (mfirst .EQ. -1 ) THEN
*       Accelerate immediately; 'no pre-acceleration' area

*        itsPreAcc = itsPreAcc + 1
*        write(*,*) "  > Not using pre-acceleration area"
        
*      ELSE
*       Find some areas BEFORE accelerating

        stopPreAccelerate = .FALSE.
        
 115    IF ( .NOT.(stopPreAccelerate) ) THEN
          itsPreAcc = itsPreAcc + 1
          write(*,*) "Using n =", itsPreAcc
          zeroStartPoint = (itsPreAcc + 1) * pi / Cy
          zeroBoundL = zeroR
          zeroBoundR = zeroStartPoint + 0.75d0 * pi / Cy

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
          write(*,*) zeroL, "and ", zeroR, "; sum: ", sum

*         STOP condition for pre-acceleration.
*         Not sure about this...
*          if ( m .EQ. (mmax - 1) ) stopPreAccelerate = .TRUE.
*          if ( itsPreAcc .GE. 5) stopPreAccelerate = .TRUE.
          
*         Stop once the real and imaginary parts have passed their final turning points
          IF (zeroR .GT. finalTP) stopPreAccelerate = .TRUE.
          
          CALL advanceM(mmax, m, mOld, leftOfMax, flip)

          GOTO 115
        ENDIF
*      ENDIF

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
*          write(*,*) "Starting with"
          itsAcceleration = itsAcceleration + 1
*         itsAcceleration = 1 means this is the first area found
*         under the acceleration regime

          n = itsPreAcc + itsAcceleration 
*          write(*,*) " Using n =", n
          zeroStartPoint = (n + 1) * pi / Cy
          zeroL = zeroR
          
          zeroBoundL = zeroStartPoint - 0.35 * pi / Cy
          zeroBoundR = zeroStartPoint + 0.15d0 * pi / Cy

          CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                        zeroStartPoint, zero)
          CALL findZeroSmallp(zero, f, df)

          zeroR = zero
          xvec(itsAcceleration + 1) = zeroR
          write(*,*) "  - Integrate between:", zeroL, zeroR

          CALL gaussq( DFintegrand, psi, zeroL, zeroR)
*         psi: area of the latest region
          wvec(itsAcceleration) = psi
          write(*,*) "  - Area between zeros is:", psi

          accMax = 40
          Wold2 = Wold
          Wold = West
          CALL accelerateNEW(xvec, wvec, itsAcceleration, 
     &                       Mmatrix, Nmatrix, West)
*          W is the best guess of the convergent integration
          write(*,*) "iteration", itsAcceleration, ":", West
          write(*,*) "  - Estimate of tail area:", West
          write(*,*) "--------------------------------"

*         Check for convergence
          write(*,*) "Ws: ", West, Wold, Wold2
          relerr = (DABS( West - Wold ) + DABS( West - Wold2 ) ) /
     &                     (DABS(West) + epsilon )
          IF (relerr .LT. aimrerr ) THEN 
            write(*,*) "  (Nice) Relerr is", relerr
            convergence = .TRUE.
          ELSE
            write(*,*) "  (Keep a-going) Relerr is", relerr
          ENDIF
        IF (itsAcceleration .EQ. 30) convergence = .TRUE.
        write(*,*) "TEMP: itsAcceleration = -30: stop"
          mOld = m
          CALL advanceM(mmax, m, mOld, leftOfMax, flip)

          GOTO 12
        ENDIF
      ELSE
        IF ( .NOT.(convergence)) THEN
        write(*,*) "  - Computing for p:", Cp
*           CALL findApproxZeros()
*           CALL integrateRegions()
*           CALL accelerate()
        ENDIF
      ENDIF
      write(*,*) "!!!!! DFsmall/big: Approx zeros can be removed !!!!!"
      
      areaA = West
      areaT = area0 + area1 + areaA
      write(*,*) "SUMMARY:"
      write(*,*) "  Area0 ", area0
      write(*,*) "  Area1 ", area1, "(", itsPreAcc, "regions)"
      write(*,*) "  AreaA ", areaA, "(", itsAcceleration, " its)"
      write(*,*) "  TOTAL ", areaT
      
      
*** WHAT TO DO with relerrr? Might have three rel errors: from initila, pre-acc, acc?
*** Take largest of the three? ADD?
*** Assume the argest relative error comes from the acceleration.
      write(*,*) "FIX rel err: |A|.relA + ... + |C|.relC/|A+B+C|"
      
*     We have the value of the integral in the CDF calculation. 
*     So now work out the CDF
      CALL findLambda(lambda)

*     The integration returns the conditional CDF for Y | Y > 0.
*     So we need to find the CDF of Y.
*     That also means adding P(Y=0) 

*     So the value returned by the integration  
     
      funvalue = -areaT/pi + 0.50d0 
      write(*,*) "FINAL AREA: The cdf value is", funvalue
      write(*,*) "DFsmallp: funvalue, exitstatus, relerr, exacti"
      write(*,*) funvalue, exitstatus, relerr, exacti

      RETURN
      END

**************************************************************************

