SUBROUTINE DFsmallp(funvalue, exitstatus, relerr, verbose) BIND(C, NAME='DFsmallp')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

! CRITICAL: IMPLICIT NONE must come directly after USE statements.
  IMPLICIT NONE

 ! --- Dummy Arguments (The variables passed into the subroutine) ---
  INTEGER(C_INT), INTENT(IN)       :: i              ! Observation index
  INTEGER(C_INT), INTENT(INOUT)    :: verbose        ! Assuming INOUT/IN for verbosity flag
  INTEGER(C_INT), INTENT(OUT)      :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: funvalue, relerr ! The final computed result and relative error

   ! --- INTERFACES: All C-bound routines called by DFbigp MUST be defined here ---
  INTERFACE
      ! 1. Function to find Kmax special point
      FUNCTION findKmaxSP(j) BIND(C, NAME='findKmaxSP')
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
          
          IMPLICIT NONE  
          REAL(KIND=C_DOUBLE)                :: findKmaxSP
          INTEGER(C_INT), INTENT(IN)  :: j
      END FUNCTION findKmaxSP

      ! 2. Subroutine to find Kmax and related indices
      SUBROUTINE findKmax(j, kmax, tmax, mmax, mfirst, startTKmax) BIND(C, NAME='findKmax')
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
        
        IMPLICIT NONE
          INTEGER(C_INT), INTENT(IN) :: j
          INTEGER(C_INT), INTENT(OUT) :: mfirst, mmax
          
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: kmax, tmax
          REAL(KIND=C_DOUBLE), INTENT(IN) :: startTKmax
      END SUBROUTINE findKmax

      ! 3. Subroutine to advance the iteration index m
      SUBROUTINE advanceM(j, m_index, mmax, mOld, leftOfMax, flip) BIND(C, NAME='advanceM')
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

          INTEGER(C_INT), INTENT(IN)      :: j, mmax
          INTEGER(C_INT), INTENT(INOUT)   :: m_index
          INTEGER(C_INT), INTENT(OUT)     :: mOld
          INTEGER(C_INT), INTENT(INOUT)   :: leftOfMax
          INTEGER(C_INT), INTENT(OUT)     :: flip
      END SUBROUTINE advanceM
      
      ! 4. Function for the integrand (used by gaussq)
      FUNCTION DFintegrand(t) BIND(C, NAME='DFintegrand')
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

          REAL(KIND=C_DOUBLE) :: DFintegrand
          REAL(KIND=C_DOUBLE), INTENT(IN) :: t
      END FUNCTION DFintegrand

      ! 5. Function for Gaussian Quadrature integration
      FUNCTION gaussq(funcd, a, b, aimrerr) BIND(C, NAME='gaussq')
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

          REAL(KIND=C_DOUBLE) :: gaussq
          INTERFACE
              FUNCTION funcd(x) BIND(C)
                USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

                REAL(KIND=C_DOUBLE) :: funcd
                REAL(KIND=C_DOUBLE), INTENT(IN) :: x
              END FUNCTION funcd
          END INTERFACE
          REAL(KIND=C_DOUBLE), INTENT(IN) :: a, b, aimrerr
      END FUNCTION gaussq

      ! 6. Subroutine for acceleration
      SUBROUTINE accelerateNEW(xvec, wvec, nzeros, Mmatrix, NMatrix, West) BIND(C, NAME='accelerateNEW')
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

          INTEGER(C_INT), INTENT(IN)      :: nzeros
          REAL(KIND=C_DOUBLE), INTENT(IN) :: xvec(200), wvec(200), Mmatrix(2, 200), Nmatrix(2, 200), West
      END SUBROUTINE accelerateNEW

      ! 7. Subroutine to find exact zeros
      SUBROUTINE findExactZeros(i, tL, tR, zeroL, zeroR) BIND(C, NAME='findExactZeros')
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
        INTEGER(C_INT), INTENT(IN)        :: i
          REAL(KIND=C_DOUBLE), INTENT(IN) :: tL, tR
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: zeroL, zeroR
      END SUBROUTINE findExactZeros

  END INTERFACE
  ! --- END INTERFACES ---

  ! --- Subroutine Outputs/Inputs ---
  REAL(KIND=C_DOUBLE) :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200), West, Wold, Wold2
  REAL(KIND=C_DOUBLE) :: zeroBoundR, zeroBoundL, zeroStartPoint

  
  ! --- Local Variables ---
  REAL(KIND=8) :: pi, area, zeroL, zeroR, zero, zeroStartPoint
  REAL(KIND=8) :: epsilon, aimrerr, omega
  INTEGER(C_INT) :: m, mmax, n, mOld
  
  INTEGER(C_INT)      :: mmax, mfirst, mOld, accMax
  INTEGER(C_INT)      :: itsacceleration, itsPreAcc, m
  INTEGER(C_INT)      :: leftOfMax, flip, convergence, stopPreAccelerate
  
  ! FIX 2: funvalue and relerr removed (already defined as dummy args)
  REAL(KIND=C_DOUBLE) :: kmax, startTKmax, tmax, aimrerr
  REAL(KIND=C_DOUBLE) :: epsilon, areaT, pi, psi, zero
  REAL(KIND=C_DOUBLE) :: zeroL, zeroR, area0, area1, areaA, sum
  REAL(KIND=8)        :: current_y, current_mu, current_phi ! Still using KIND=8 for internal module array access
  REAL(KIND=C_DOUBLE) :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200), West, Wold, Wold2
  REAL(KIND=C_DOUBLE) :: zeroBoundR, zeroBoundL, zeroStartPoint

  ! --- External Function/Subroutine Declarations (CRITICAL: only EXTERNAL needed) ---
  ! DFsmallp only calls DFintegrand, gaussq, and findExactZeros directly.
  REAL(KIND=8) DFintegrand
  
  ! --- Initialization ---
  CpSmall = .TRUE.
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
 
  pi = DACOS(1.0_8)
  exitstatus = 0
  relerr = 0.0_8 
  epsilon = 1.0_8 / 10.0_8**16 ! 1.0d-16
  aimrerr = 1.0_8 / 10.0_8**14 ! 1.0d-14
  
     IF (verbose .EQ. 1) WRITE(*,*) " FOR 1 < p < 2"
      ! Initialise the M and N matrices, x and w vectors
      DO i = 1, 200
        Mmatrix(1, i) = 0.0d0
        Mmatrix(2, i) = 0.0d0
        Nmatrix(1, i) = 0.0d0
        Nmatrix(2, i) = 0.0d0
        xvec(i) = 0.d000
        wvec(i) = 0.0d00
      ENDDO
      
      ! Set up
      pi = 4.0d0 * DATAN(1.0d0)
      exitstatus = 0
      relerr = 1.0d00
      epsilon = 1.0d-16
      aimrerr = 1.0d-14
      convergence = 0


      ! FIND kmax, tmax, mmax
      IF (Cy(i) .GE. Cmu(i)) THEN
        ! ************** y > MU   **************
        ! Im k(t) heads down immediately

        IF (verbose .EQ. 1) WRITE(*,*) "** y >= mu"
        kmax = 0.0d00
        tmax = 0.0d00
        mmax = 0
        mfirst = -1
        mOld = 0
      
        zeroStartPoint = pi / Cy
        leftOfMax = .FALSE.
      ELSE
        ! ************** y < MU   **************
        ! HARDER!         
        IF (verbose .EQ. 1) WRITE(*,*) "** y < mu"
        
        WRITE(*,*) "About to find kmax"
        startTKMax = findKmaxSP(i)
        IF (verbose .EQ. 1) WRITE(*,*) "Find kmax, start at: ", StartTKmax
        
        ! Sometimes, important to spend some getting a good starting point and bounds.  
        CALL findKmaxSPbounds(startTKmax, kmaxL, kmaxR)

        startTKmax =  (kmaxL + kmaxR) / 2.0d0
        
        leftOfMax = 1
        IF ( mmax .EQ. 0) THEN
          mfirst = 0
          mOld = 0
          zeroStartPoint = tmax + pi/Cy(i)
          leftOfMax = 0
        ELSE
          mfirst = 1
          mOld = 0
          zeroStartPoint = pi / (Cmu(i) - Cy(i))
          mOld = m

          CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

        END IF
        IF (verbose .EQ. 1) WRITE(*,*) "Find kmax, start at: ", StartTKmax
        CALL findKmax(i,kmax, tmax, mmax, mfirst, startTKmax, kmaxL, kmaxR)
  
        IF (verbose .EQ. 1) THEN
          WRITE(*,*) "** Found(b): kmax =", kmax
          WRITE(*,*) "             tmax =", tmax
          WRITE(*,*) "             mmax =", mmax
        END IF
        
        leftOfMax = 1
        IF ( mmax .EQ. 0) THEN
          mfirst = 0
          mOld = 0
          zeroStartPoint = tmax + pi/Cy(i)
          leftOfMax = 0
        ELSE
          mfirst = 1
          mOld = 0
          zeroStartPoint = pi / (Cmu(i) - Cy(i))
          mOld = m

          CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
        END IF
      
      END IF

      WRITE(*,*) "--- (Deal with returned errors, non-convergence)"

      ! INTEGRATION
      ! There are three integration regions:
      !   1. The *initial* area, which is not between zeros of Im{k(t)}: area0
      !   2. The initial area *before* Sidi acceleration is invoked: area1
      !      (For instance, wait until Im{k(t)} is on the downturn.)
      !   3. The area thereafter, upon which Sidi acceleration is
      !      applied; the area returned by acceleration is areaA

      area0 = 0.0d00
      area1 = 0.0d00
      areaA = 0.0d00
      m = mfirst
      
      ! Find the final turning point of Im/Re k, and start accelerating thereafter      
      WRITE(*,*) "IS this about TPs correct??"
      finalTP = 0.0d00
      IF ( Cy .LT. Cmu ) THEN
        CALL findAccelStart(finalTP)
      END IF
      
      zeroStartPoint = pi / Cy
      ! TRY A NEW ONE!
      front = Cmu ** (1.0d0 - Cp) / ( Cphi * (1.0d0 - Cp))
      zeroStartPoint = front * DTAN( pi * ( 1.0d0 - Cp) / Cp )
      WRITE(*,*) "zeroStartPoint", zeroStartPoint

      ! 1. INTEGRATE FIRST REGION: area0
      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "*******************************" 
        WRITE(*,*) "1. INTEGRATE: the INITIAL region"
      END IF 
      zeroBoundL = tmax
      zeroBoundR = zeroStartPoint + 0.25d0 * pi / Cy(i)
      IF (verbose .EQ. 1) WRITE(*,*) " Bounds zero; ", zeroBoundL, zeroBoundR
      WRITE(*,*) "   m = ", m

      ! Now find the right-side zero
      CALL findExactZeros(zeroBoundL, zeroBoundR, zeroStartPoint, zero)

      zeroL =  0.0d00
      zeroR = zero
      area0 = gaussq(area0, zeroL, zeroR, aimrerr)

      IF (verbose .EQ. 1) WRITE(*,*) "  - Initial area:", area0
      IF (verbose .EQ. 1) WRITE(*,*) "    between:", zeroL, zeroR


      ! 2. INTEGRATE: the PRE-ACCELERATION regions: area1
      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "*******************************" 
        WRITE(*,*) "2. INTEGRATE: Pre-acceleration"
      END IF 

      ! When p > 2, things seem well-behaved most of the time, so 
      ! we declare  area1  to be up to m = mmax - 1 (i.e., just
      ! after the downturn)

      itsPreAcc = 0
      area1 = 0.0d00
      CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

      ! IF (mfirst .EQ. -1 ) THEN
      !  Accelerate immediately; 'no pre-acceleration' area

      ! itsPreAcc = itsPreAcc + 1
      ! WRITE(*,*) "  > Not using pre-acceleration area"
        
      ! ELSE
      ! Find some areas BEFORE accelerating

        stopPreAccelerate = .FALSE.
        
 115    IF ( .NOT.(stopPreAccelerate) ) THEN
          itsPreAcc = itsPreAcc + 1
          zeroL = zeroR

          zeroStartPoint = (itsPreAcc + 1) * pi / Cy(i)
          ! WRITE(*,*)" StartPT:", zeroStartPoint
          zeroBoundL = zeroR
          zeroBoundR = zeroStartPoint + 0.75d0 * pi / Cy(i)

          CALL findExactZeros(zeroBoundL, zeroBoundR, zeroStartPoint, zero)

          zeroR = zero

          CALL gaussq( DFintegrand, sum, zeroL, zeroR)
          area1 = area1 + sum
          
          IF (verbose .EQ. 1) THEN
            WRITE(*,*) "--- Integrate (m = ", m, ") between " 
            WRITE(*,*) zeroL, "and ", zeroR, "; sum: ", sum
          END IF

          ! STOP condition for pre-acceleration.
          ! Not sure about this...
          ! if ( m .EQ. (mmax - 1) ) stopPreAccelerate = .TRUE.
          ! if ( itsPreAcc .GE. 5) stopPreAccelerate = .TRUE.
          
          ! Stop once the real and imaginary parts have passed their final turning points
          IF (zeroR .GT. finalTP) stopPreAccelerate = 1
          
          CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

          GOTO 115
        END IF
        ! END IF

      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "Finished pre-acc; areas"
        WRITE(*,*) "  Area1 ", area1
      END IF


      ! 3. INTEGRATE: the ACCELERATION regions: areaA
      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "*******************************" 
        WRITE(*,*) "3. INTEGRATE: the ACCELERATION"
      END IF
      
      Wold = 0.0d00
      Wold2 = 1.0d00


      itsAcceleration = 0
      areaA = 0.0d00
      convergence = 0

      xvec(1) = zeroR
      ! This will be the very first, left-most value of t used, the left-most
      ! value of  t  used in the acceleration (the previous regions *right* value) 

 12   IF ( convergence .EQ. 0) THEN
        itsAcceleration = itsAcceleration + 1
        ! itsAcceleration = 1 means this is the first area found
        ! under the acceleration regime

        n = itsPreAcc + itsAcceleration 
        zeroStartPoint = (n + 1) * pi / Cy(i)
        zeroL = zeroR
          
        zeroBoundL = zeroStartPoint - 0.35 * pi / Cy(i)
        zeroBoundR = zeroStartPoint + 0.15d0 * pi / Cy(i)

        CALL findExactZeros(zeroBoundL, zeroBoundR, zeroStartPoint, zero)
        CALL findZeroSmallp(zero, f, df)

        zeroR = zero
        xvec(itsAcceleration + 1) = zeroR
        IF (verbose .EQ. 1) WRITE(*,*) "  - Integrate between:", zeroL, zeroR

        CALL gaussq( DFintegrand, psi, zeroL, zeroR)
        ! psi: area of the latest region
        wvec(itsAcceleration) = psi
        IF (verbose .EQ. 1) WRITE(*,*) "  - Area between zeros is:", psi

        accMax = 100
        Wold2 = Wold
        Wold = West
        CALL accelerateNEW(xvec, wvec, itsAcceleration, Mmatrix, Nmatrix, West)
        ! W is the best guess of the convergent integration
         if (verbose .EQ. 1) WRITE(*,*) "  - Tail estimate:", West

      ! Check for convergence
        relerr = (DABS( West - Wold ) + DABS( West - Wold2 ) ) / (DABS(West) + epsilon )
        IF (relerr .LT. aimrerr ) THEN 
          convergence = 1
        END IF

        IF (itsAcceleration .EQ. accMax) THEN
          convergence = 1
          WRITE(*,*) "No convergence of acceleration."
        END IF

        mOld = m
        CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

        GOTO 12
      END IF
      WRITE(*,*) "!!!!! DFsmall/big: Approx zeros can be removed !!!!!"
      
      areaA = West
      areaT = area0 + area1 + areaA
      
      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "SUMMARY:"
        WRITE(*,*) "  * Area0 ", area0
        WRITE(*,*) "  * Area1 ", area1, "(", itsPreAcc, "regions)"
        WRITE(*,*) "  * AreaA ", areaA, "(", itsAcceleration, " its)"
        WRITE(*,*) "  TOTAL ", areaT
      END IF
      
      ! WHAT TO DO with relerrr? Might have three rel errors: from initila, pre-acc, acc?
      ! Take largest of the three? ADD?
      ! Assume the argest relative error comes from the acceleration.
      WRITE(*,*) "FIX rel err: |A|.relA + ... + |C|.relC/|A+B+C|"
      
      ! We have the value of the integral in the CDF calculation. 
      ! So now work out the CDF
      CALL findLambda(lambda)

      ! The integration returns the conditional CDF for Y | Y > 0.
      ! So we need to find the CDF of Y.
      ! That also means adding P(Y=0) 

      ! So the value returned by the integration  
     
      funvalue = -areaT/pi + 0.50d0 
      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "FINAL AREA: The cdf value is", funvalue
        WRITE(*,*) "DFsmallp: funvalue, exitstatus, relerr"
        WRITE(*,*) funvalue, exitstatus, relerr
      END IF
      
END SUBROUTINE DFsmallp
