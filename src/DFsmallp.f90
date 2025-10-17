SUBROUTINE DFsmallp(i, funvalue, exitstatus, relerr, verbose) 
  USE tweedie_params_mod

  IMPLICIT NONE

 ! --- Dummy Arguments, variables passed into the subroutine
  INTEGER, INTENT(IN)                      :: i               ! Observation index
  REAL(KIND=8), DIMENSION(*), INTENT(OUT)  :: funvalue        ! The final computed result and relative error
  INTEGER, INTENT(OUT)                     :: exitstatus      ! Output status
  REAL(KIND=8), INTENT(OUT)                :: relerr          ! The final computed result and relative error
  INTEGER, INTENT(INOUT)                   :: verbose         ! Assuming INOUT/IN for verbosity flag

   ! --- INTERFACES: All C-bound routines called by DFsmallp:
  INTERFACE
    FUNCTION findKmaxSP(j)
      IMPLICIT NONE  
      REAL(KIND=8)          :: findKmaxSP
      INTEGER, INTENT(IN)   :: j
    END FUNCTION findKmaxSP


    SUBROUTINE findKmax(j, kmax, tmax, mmax, mfirst, startTKmax)
      IMPLICIT NONE
      INTEGER, INTENT(IN)         :: j
      INTEGER, INTENT(OUT)        :: mfirst, mmax
      REAL(KIND=8), INTENT(OUT)   :: kmax, tmax
      REAL(KIND=8), INTENT(IN)    :: startTKmax
    END SUBROUTINE findKmax


    SUBROUTINE advanceM(j, m_index, mmax, mOld, leftOfMax, flip) 
      INTEGER, INTENT(IN)       :: j, mmax
      INTEGER, INTENT(INOUT)    :: m_index
      INTEGER, INTENT(OUT)      :: mOld
      INTEGER, INTENT(INOUT)    :: leftOfMax
      INTEGER, INTENT(OUT)      :: flip
    END SUBROUTINE advanceM
      

    SUBROUTINE DFgaussq(i, a, b, area)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)          :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: a, b
    END SUBROUTINE DFgaussq


    SUBROUTINE acceleratenew(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)      :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(IN) :: xvec(200), wvec(200), Mmatrix(2, 200), Nmatrix(2, 200), West
    END SUBROUTINE acceleratenew


    SUBROUTINE findExactZeros(i, m, tL, tR, zeroL, zeroR)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i, m
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: tL, tR
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: zeroL, zeroR
    END SUBROUTINE findExactZeros


    SUBROUTINE findZeroSmallp(i, t, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
    END SUBROUTINE findZeroSmallp


    SUBROUTINE findLambda(i, lambda)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: lambda
    END SUBROUTINE findLambda
  END INTERFACE
  ! --- END INTERFACES ---

  ! --- Subroutine Outputs/Inputs ---
  REAL(KIND=C_DOUBLE) :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200), West, Wold, Wold2
  REAL(KIND=C_DOUBLE) :: zeroBoundR, zeroBoundL, zeroStartPoint

  
  ! --- Local Variables ---
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR, zero, aimrerr
  INTEGER               :: m, n, mOld, mmax, mfirst, accMax, j
  INTEGER               :: leftOfMax, flip, convergence, stopPreAccelerate

  REAL(KIND=8)          :: pi
  INTEGER               :: itsacceleration, itsPreAcc

  ! FIX 2: funvalue and relerr removed (already defined as dummy args)
  REAL(KIND=C_DOUBLE)  :: kmax, startTKmax, tmax, df, f, finalTP, front, kmaxL, kmaxR, lambda
  REAL(KIND=C_DOUBLE)  :: areaT, area0, area1, areaA, sum, psi, epsilon
  REAL(KIND=C_DOUBLE)  :: current_y, current_mu, current_phi ! Still using KIND=8 for internal module array access

  ! --- Initialization ---

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
 
  exitstatus = 0
  relerr = 0.0_8 
  epsilon = 1.0_8 / 10.0_8**16 ! 1.0d-16
  aimrerr = 1.0_8 / 10.0_8**14 ! 1.0d-14

  IF (verbose .EQ. 1) WRITE(*,*) " FOR 1 < p < 2"
    ! Initialise the M and N matrices, x and w vectors
    DO j = 1, 200
      Mmatrix(1, j) = 0.0d0
      Mmatrix(2, j) = 0.0d0
      Nmatrix(1, j) = 0.0d0
      Nmatrix(2, j) = 0.0d0
      xvec(j) = 0.d000
      wvec(j) = 0.0d00
    END DO
    
    ! Set up
    pi = 4.0D0 * DATAN(1.0D0)
    exitstatus = 0
    relerr = 1.0d00
    epsilon = 1.0d-16
    aimrerr = 1.0d-14
    convergence = 0
    mOld = 0
    m = 0

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
    
      zeroStartPoint = pi / Cy(i)
      leftOfMax = 0
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
      
        IF (verbose .EQ. 1) WRITE(*,*) "Find kmax, start at: ", StartTKmax
        CALL findKmax(i, kmax, tmax, mmax, mfirst, startTKmax)
          !  WAS:            CALL findKmax(i, kmax, tmax, mmax, mfirst, startTKmax, kmaxL, kmaxR)

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
      IF ( Cy(i) .LT. Cmu(i) ) THEN
        CALL findAccelStart(finalTP)
      END IF
    
      zeroStartPoint = pi / Cy(i)
      ! TRY A NEW ONE!
      front = Cmu(i) ** (1.0d0 - Cp) / ( Cphi(i) * (1.0d0 - Cp))
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
      CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)

      zeroL =  0.0d00
      zeroR = zero
      CALL DFgaussq(i, zeroL, zeroR, area0)

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

      stopPreAccelerate = 0
        
      DO WHILE ( stopPreAccelerate .EQ. 0 )
        itsPreAcc = itsPreAcc + 1
        zeroL = zeroR

        zeroStartPoint = (itsPreAcc + 1) * pi / Cy(i)
        ! WRITE(*,*)" StartPT:", zeroStartPoint
        zeroBoundL = zeroR
        zeroBoundR = zeroStartPoint + 0.75d0 * pi / Cy(i)

        CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)

        zeroR = zero

        CALL DFgaussq(i, zeroL, zeroR, sum)
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

      END DO
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

      DO WHILE ( convergence .EQ. 0)
        itsAcceleration = itsAcceleration + 1
        ! itsAcceleration = 1 means this is the first area found
        ! under the acceleration regime

        n = itsPreAcc + itsAcceleration 
        zeroStartPoint = (n + 1) * pi / Cy(i)
        zeroL = zeroR
          
        zeroBoundL = zeroStartPoint - 0.35 * pi / Cy(i)
        zeroBoundR = zeroStartPoint + 0.15d0 * pi / Cy(i)

    WRITE(*,*) "ABOU TO CALL IN DFbigf"
        CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)
    WRITE(*,*) "BACK IN DFbigf; about to call findZeroSmallp"
        CALL findZeroSmallp(i, zero, f, df)
    WRITE(*,*) "OUT FROM findZeroSmallp"

        zeroR = zero
        xvec(itsAcceleration + 1) = zeroR
        IF (verbose .EQ. 1) WRITE(*,*) "  - Integrate between:", zeroL, zeroR

        CALL DFgaussq(i, zeroL, zeroR, psi)
        ! psi: area of the latest region
        wvec(itsAcceleration) = psi
        IF (verbose .EQ. 1) WRITE(*,*) "  - Area between zeros is:", psi

        accMax = 100
        Wold2 = Wold
        Wold = West
        CALL acceleratenew(xvec, wvec, itsAcceleration, Mmatrix, Nmatrix, West)
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

      END DO
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
      CALL findLambda(i, lambda)

      ! The integration returns the conditional CDF for Y | Y > 0.
      ! So we need to find the CDF of Y.
      ! That also means adding P(Y=0) 

      ! So the value returned by the integration  
     
      funvalue(i) = -areaT/pi + 0.50d0 
      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "FINAL AREA: The cdf value is", funvalue(i)
        WRITE(*,*) "DFsmallp: funvalue, exitstatus, relerr"
        WRITE(*,*) funvalue(i), exitstatus, relerr
      END IF
    
END SUBROUTINE DFsmallp
