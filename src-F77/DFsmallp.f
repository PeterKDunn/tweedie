
      SUBROUTINE DFsmallp(funvalue, exitstatus, relerr, verbose)

*     Calculates the DF of the log-likelihood function of a
*     Poisson-gamma distribution by inverting the MGF: 1 < p < 2
*
*     OUT:  funvalue, exitstatus, relerr, its

      IMPLICIT NONE
      DOUBLE PRECISION funvalue, pi, zero, zeroL, zeroR, sum
      DOUBLE PRECISION aimrerr, relerr, tmax, kmax, f, df
      DOUBLE PRECISION kmaxL, kmaxR
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      DOUBLE PRECISION findKmaxSP, startTKMax, front
      DOUBLE PRECISION zeroStartPoint, area0, area1, areaA
      DOUBLE PRECISION zeroBoundL, zeroBoundR, Integrands, psi
      DOUBLE PRECISION Wold, Wold2, areaT, epsilon
      DOUBLE PRECISION West, xvec(300), wvec(300), lambda
      DOUBLE PRECISION Mmatrix(2, 200), Nmatrix(2, 200)
      DOUBLE PRECISION finalTP
      INTEGER mfirst, m, mOld, exitstatus, mmax, n, i
      INTEGER itsPreAcc, accMax, itsAcceleration
      LOGICAL  convergence, flip, leftOfMax
      LOGICAL pSmall, stopPreAccelerate, verbose
      EXTERNAL findKmaxSP, Integrands
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

      IF (verbose) write(*,*) " FOR 1 < p < 2"
      
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


*     FIND kmax, tmax, mmax
      IF (Cy. GE. Cmu) THEN
*     ************** y > MU   **************
*     Im k(t) heads down immediately

        IF (verbose) write(*,*) "** y >= mu"
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
        IF (verbose) write(*,*) "** y < mu"
        
        write(*,*) "About to find kmax"
        startTKMax = findKmaxSP()
        IF (verbose) write(*,*) "Find kmax, start at: ", StartTKmax
        
*       Sometimes, important to spend some getting a good starting point and bounds.  
        CALL improveKmaxSPBounds(startTKmax, kmaxL, kmaxR)

        startTKmax =  (kmaxL + kmaxR) / 2.0d0
        
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
        IF (verbose) write(*,*) "Find kmax, start at: ", StartTKmax
        CALL findKmax(kmax, tmax, mmax, mfirst, startTKmax,
     &                kmaxL, kmaxR)
  
        IF (verbose) THEN
          write(*,*) "** Found(b): kmax =", kmax
          write(*,*) "             tmax =", tmax
          write(*,*) "             mmax =", mmax
        ENDIF
        
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
      m = mfirst
      
*     Find the final turning point of Im/Re k, and start accelerating thereafter      
      write(*,*) "IS this about TPs correct??"
      finalTP = 0.0d00
      IF ( Cy .LT. Cmu ) THEN
        CALL findAccelStart(finalTP)
      ENDIF
      
      zeroStartPoint = pi / Cy
* TRY A NEW ONE!
      front = Cmu ** (1.0d0 - Cp) / ( Cphi * (1.0d0 - Cp))
      zeroStartPoint = front * DTAN( pi * ( 1.0d0 - Cp) / Cp )
      write(*,*) "zeroStartPoint", zeroStartPoint

*     1. INTEGRATE FIRST REGION: area0
      IF (verbose) THEN
        write(*,*) "*******************************" 
        write(*,*) "1. INTEGRATE: the INITIAL region"
      ENDIF 
      zeroBoundL = tmax
      zeroBoundR = zeroStartPoint + 0.25d0 * pi / Cy
      IF (verbose) write(*,*) " Bounds zero; ", zeroBoundL, zeroBoundR
      write(*,*) "   m = ", m

*     Now find the right-side zero
      CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                    zeroStartPoint, zero)

      zeroL =  0.0d00
      zeroR = zero

      CALL gaussq( Integrands, area0, zeroL, zeroR)
      IF (verbose) write(*,*) "  - Initial area:", area0
      IF (verbose) write(*,*) "    between:", zeroL, zeroR


*     2. INTEGRATE: the PRE-ACCELERATION regions: area1
      IF (verbose) THEN
        write(*,*) "*******************************" 
        write(*,*) "2. INTEGRATE: Pre-acceleration"
      ENDIF 

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
          zeroL = zeroR

          zeroStartPoint = (itsPreAcc + 1) * pi / Cy
*          write(*,*)" StartPT:", zeroStartPoint
          zeroBoundL = zeroR
          zeroBoundR = zeroStartPoint + 0.75d0 * pi / Cy

          CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                        zeroStartPoint, zero)

          zeroR = zero

          CALL gaussq( Integrands, sum, zeroL, zeroR)
          area1 = area1 + sum
          
          IF (verbose) THEN
            write(*,*) "--- Integrate (m = ", m, ") between " 
            write(*,*) zeroL, "and ", zeroR, "; sum: ", sum
          ENDIF

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

      IF (verbose) THEN
        write(*,*) "Finished pre-acc; areas"
        write(*,*) "  Area1 ", area1
      ENDIF


*     3. INTEGRATE: the ACCELERATION regions: areaA
      IF (verbose) THEN
        write(*,*) "*******************************" 
        write(*,*) "3. INTEGRATE: the ACCELERATION"
      ENDIF
      
      Wold = 0.0d00
      Wold2 = 1.0d00


      itsAcceleration = 0
      areaA = 0.0d00
      convergence = .FALSE.

      xvec(1) = zeroR
*     This will be the very first, left-most value of t used, the left-most
*     value of  t  used in the acceleration (the previous regions *right* value) 

 12   IF ( .NOT.(convergence)) THEN
        itsAcceleration = itsAcceleration + 1
*       itsAcceleration = 1 means this is the first area found
*       under the acceleration regime

        n = itsPreAcc + itsAcceleration 
        zeroStartPoint = (n + 1) * pi / Cy
        zeroL = zeroR
          
        zeroBoundL = zeroStartPoint - 0.35 * pi / Cy
        zeroBoundR = zeroStartPoint + 0.15d0 * pi / Cy

        CALL findExactZeros(zeroBoundL, zeroBoundR, 
     &                      zeroStartPoint, zero)
        CALL findZeroSmallp(zero, f, df)

        zeroR = zero
        xvec(itsAcceleration + 1) = zeroR
        IF (verbose) write(*,*) "  - Integrate between:", zeroL, zeroR

        CALL gaussq( Integrands, psi, zeroL, zeroR)
*       psi: area of the latest region
        wvec(itsAcceleration) = psi
        IF (verbose) write(*,*) "  - Area between zeros is:", psi

        accMax = 100
        Wold2 = Wold
        Wold = West
        CALL accelerateNEW(xvec, wvec, itsAcceleration, 
     &                     Mmatrix, Nmatrix, West)
*        W is the best guess of the convergent integration
         if (verbose) write(*,*) "  - Tail estimate:", West

*       Check for convergence
        relerr = (DABS( West - Wold ) + DABS( West - Wold2 ) ) /
     &                   (DABS(West) + epsilon )
        IF (relerr .LT. aimrerr ) THEN 
          convergence = .TRUE.
        ENDIF

        IF (itsAcceleration .EQ. accMax) THEN
          convergence = .TRUE.
          write(*,*) "No convergence of acceleration."
        ENDIF

        mOld = m
        CALL advanceM(mmax, m, mOld, leftOfMax, flip)

        GOTO 12
      ENDIF
      write(*,*) "!!!!! DFsmall/big: Approx zeros can be removed !!!!!"
      
      areaA = West
      areaT = area0 + area1 + areaA
      
      IF (verbose) THEN
        write(*,*) "SUMMARY:"
        write(*,*) "  * Area0 ", area0
        write(*,*) "  * Area1 ", area1, "(", itsPreAcc, "regions)"
        write(*,*) "  * AreaA ", areaA, "(", itsAcceleration, " its)"
        write(*,*) "  TOTAL ", areaT
      ENDIF
      
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
      IF (verbose) THEN
        write(*,*) "FINAL AREA: The cdf value is", funvalue
        write(*,*) "DFsmallp: funvalue, exitstatus, relerr"
        write(*,*) funvalue, exitstatus, relerr
      ENDIF
      
      RETURN
      END

**************************************************************************

