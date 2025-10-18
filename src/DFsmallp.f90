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


    SUBROUTINE findKmaxSPbounds(i, startx, xL, xR)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: startx
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: xL, xR
    END SUBROUTINE findKmaxSPbounds


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


    SUBROUTINE findAccelStart(i, tRight)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: tRight
    END SUBROUTINE findAccelStart


    SUBROUTINE findLambda(i, lambda)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: lambda
    END SUBROUTINE findLambda
  END INTERFACE
  ! --- END INTERFACES ---

  ! --- Subroutine Outputs/Inputs ---
  REAL(KIND=C_DOUBLE)             :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200)
  REAL(KIND=C_DOUBLE), VOLATILE   :: West, Wold, Wold2
  REAL(KIND=C_DOUBLE)             :: zeroBoundR, zeroBoundL, zeroStartPoint

  
  ! --- Local Variables ---
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR, zero, aimrerr, epsilon
  INTEGER               :: m, n, mOld, mmax, mfirst, accMax, j, minAccRegions
  INTEGER               :: leftOfMax, flip, convergence, stopPreAccelerate

  REAL(KIND=8)          :: pi
  INTEGER               :: itsacceleration, itsPreAcc

  ! FIX 2: funvalue and relerr removed (already defined as dummy args)
  REAL(KIND=C_DOUBLE)  :: kmax, startTKmax, tmax, df, f, finalTP, front, kmaxL, kmaxR, lambda
  REAL(KIND=C_DOUBLE)  :: areaT, area0, area1, areaA, sum, psi
  REAL(KIND=C_DOUBLE)  :: current_y, current_mu, current_phi ! Still using KIND=8 for internal module array access

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
 
  ! --- Initialization ---
  pi = 4.0D0 * DATAN(1.0D0)
  aimrerr = 1.0D-12
  mOld = 0
  m = 0
  exitstatus = 0
  relerr = 1.0_8
  convergence = 0
  epsilon = 1.0d-12


  IF (verbose .EQ. 1) WRITE(*,*) " FOR 1 < p < 2"
    ! Initialise the M and N matrices, x and w vectors
!    DO j = 1, 200
!      Mmatrix(1, j) = 0.0d0
!      Mmatrix(2, j) = 0.0d0
!      Nmatrix(1, j) = 0.0d0
!      Nmatrix(2, j) = 0.0d0
!      xvec(j) = 0.d000
!      wvec(j) = 0.0d00
!    END DO

    ! FIND kmax, tmax, mmax
  IF (current_y .GE. current_mu) THEN
    ! ************** y > MU   **************
    ! Im k(t) heads down immediately

    IF (verbose .EQ. 1) WRITE(*,*) "** y >= mu"
    kmax = 0.0d00
    tmax = 0.0d00
    mmax = 0
    mfirst = -1
    mOld = 0

    IF (verbose .EQ. 1) WRITE(*,*) "** Im k(t) heads down immediately"
        
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
  WRITE(*,*) "About to find SP bounds"
    CALL findKmaxSPbounds(i, startTKmax, kmaxL, kmaxR)
  WRITE(*,*) "FOUND SP bounds"
    startTKmax =  (kmaxL + kmaxR) / 2.0d0
      
    IF (verbose .EQ. 1) WRITE(*,*) "Find kmax, start at: ", StartTKmax
    CALL findKmax(i, kmax, tmax, mmax, mfirst, startTKmax)
    !  WAS:            CALL findKmax(i, kmax, tmax, mmax, mfirst, startTKmax, kmaxL, kmaxR)
    WRITE(*,*) "FINDKMAX: kmax, tmax, mmax: ", kmax, tmax, mmax

    leftOfMax = 1
    IF ( mmax .EQ. 0) THEN
      mfirst = 0
      mOld = 0
      zeroStartPoint = tmax + pi/current_y
      leftOfMax = 0
    ELSE
      mfirst = 1
      mOld = 0
      zeroStartPoint = pi / (current_mu - current_y)
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
  
  WRITE(*,*) "** Found(b): kmax =", kmax
  WRITE(*,*) "             tmax =", tmax
  WRITE(*,*) "             mmax =", mmax
  WRITE(*,*) "--- (Deal with returned errors, non-convergence)"

  ! INTEGRATION
  ! There are three integration regions:
  !   1. The *initial* area, which is not between zeros of Im{k(t)}: area0
  !   2. The initial area *before* Sidi acceleration is invoked: area1
  !      (For instance, wait until Im{k(t)} is on the downturn.)
  !   3. The area thereafter, upon which Sidi acceleration is
  !      applied; the area returned by acceleration is areaA

  ! --- Integration initialization ---
  area0 = 0.0d00
  area1 = 0.0d00
  areaA = 0.0d00
  m = mfirst
  ! Find the final turning point of Im/Re k, and start accelerating thereafter      
  WRITE(*,*) "IS this about TPs correct??"
  finalTP = 0.0d00
  IF ( current_y .LT. current_mu ) THEN
    ! Find where the oscillations settle down, so we can use acceleration
    CALL findAccelStart(i, finalTP)
    WRITE(*,*) "Start acc at: ", finalTP
  END IF

  zeroStartPoint = pi / current_y
  ! TRY A NEW ONE!
  front = current_mu ** (1.0d0 - Cp) / ( current_phi * (1.0d0 - Cp))
  zeroStartPoint = front * DTAN( pi * ( 1.0d0 - Cp) / Cp )
  write(*,*) "zeroStartPoint", zeroStartPoint
  

  ! 1. INTEGRATE FIRST REGION: area0
  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "*******************************" 
    WRITE(*,*) "1. INTEGRATE: the INITIAL region"
  END IF 
  zeroBoundL = tmax
  zeroBoundR = zeroStartPoint + 3.0d0 * pi / current_y
  IF (verbose .EQ. 1) WRITE(*,*) " Bounds zero; ", zeroBoundL, zeroBoundR

  CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)
  zeroL = 0.0d00
  zeroR = zero

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

    zeroStartPoint = (itsPreAcc + 1) * pi / current_y
    ! WRITE(*,*)" StartPT:", zeroStartPoint
    zeroBoundL = zeroR
    zeroBoundR = zeroStartPoint + 0.75d0 * pi / current_y

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

  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "Finished pre-acc; areas"
    WRITE(*,*) "  Area1 ", area1
  END IF


  ! 3. INTEGRATE: the ACCELERATION regions: areaA
  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "*******************************" 
    WRITE(*,*) "3. INTEGRATE: the ACCELERATION"
  END IF
      
 ! Initialisation
  West = 3.0_8
  Wold = 2.0_8
  Wold2 = 1.0_8
  itsAcceleration = 0
  areaA = 0.0_8
  convergence = 0
  accMax = 40           ! Maximum number of regions in accelerationl arbitrary
  minAccRegions = 3     ! Minimum number of acceleration regions to use
  

  xvec(1) = zeroR
  ! This will be the very first, left-most value of t used, the left-most
  ! value of  t  used in the acceleration (the previous regions *right* value) 

      DO WHILE ( convergence .EQ. 0)
    IF (verbose .EQ. 1) WRITE(*,*) "  --------------------- "
    IF (verbose .EQ. 1) WRITE(*,*) "  Next tail region"

        itsAcceleration = itsAcceleration + 1

!!! SUFELY WE ARE ON TEH RIGHT OF TMAX!!!!
    IF (leftOfMax .EQ. 1) THEN
      zeroStartPoint = zeroR
      zeroL = zeroR
      zeroR = zeroR * 20.0_8
    ELSE
      IF (flip .EQ. 1) THEN
        ! FLIPPING to other side of tmax
        zeroStartPoint = tmax + (tmax - zero)
        ! That is, start of the other side of tmax
        zeroL = zero
        zeroR = zeroStartPoint * 20.0_8
      ELSE
        zeroStartPoint = zeroR
        zeroL = zeroR
        zeroR = zeroR * 10.0_8
      END IF
    END IF

!        n = itsPreAcc + itsAcceleration 
!        zeroStartPoint = (n + 1) * pi / current_y
!        zeroL = zeroR

!        zeroBoundL = zeroStartPoint - 0.35d0 * pi / current_y
!        zeroBoundR = zeroStartPoint + 0.15d0 * pi / current_y
    WRITE(*,*) "  -> updated bounds:", zeroBoundL, zeroBoundR
    WRITE(*,*) "     updated start pt:", zeroStartPoint
        CALL findExactZeros(i, m, zeroL, zeroR, zeroStartPoint, zero)
    WRITE(*,*) " ::::::: zero", zero
! WHAT IS THIS NEXT LINE FOR?????????????????????????????????????????????????
        CALL findZeroSmallp(i, zero, f, df)
    WRITE(*,*) " ::::::: zero", zero

        zeroR = zero
        xvec(itsAcceleration + 1) = zeroR
        IF (verbose .EQ. 1) WRITE(*,*) "  - Integrate between:", zeroL, zeroR
        IF (verbose .EQ. 1) WRITE(*,*) "    : where m = ", m

        CALL DFgaussq(i, zeroL, zeroR, psi)
        ! psi: area of the latest region
        wvec(itsAcceleration) = psi
        IF (verbose .EQ. 1) WRITE(*,*) "  - Area between zeros is:", psi

    Wold2 = Wold
    Wold = West
    CALL acceleratenew(xvec, wvec, itsAcceleration, Mmatrix, Nmatrix, West)

        ! W is the best guess of the convergent integration
         if (verbose .EQ. 1) WRITE(*,*) "  - Tail estimate:", West

      ! Check for convergence
        relerr = (DABS( West - Wold ) + DABS( West - Wold2 ) ) / (DABS(West) + epsilon )
    ! Declare convergence of we have sufficient regions, and relerr estimate is small
    IF ( (itsAcceleration .GE. minAccRegions) .AND. &
         (relerr .LT. aimrerr) ) THEN
      IF (verbose .EQ. 1) WRITE(*,*) "  Relerr is", relerr
      convergence = 1
    END IF
    WRITE(*,*) "CONVERGENCE", convergence
    mOld = m


        IF (itsAcceleration .EQ. accMax) THEN
          convergence = 1
          WRITE(*,*) "No convergence of acceleration."
        END IF

        mOld = m
        CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

      END DO

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
        WRITE(*,*) "EXP(-lanbda) = ", DEXP(-lambda)
        WRITE(*,*) "DFsmallp: funvalue, exitstatus, relerr"
        WRITE(*,*) funvalue(i), exitstatus, relerr
      END IF
    
END SUBROUTINE DFsmallp
