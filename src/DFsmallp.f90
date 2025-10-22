SUBROUTINE DFsmallp(i, funvalueI, exitstatus, relerr, verbose, count_Integration_Regions)  
  USE tweedie_params_mod

  IMPLICIT NONE

 ! --- Dummy Arguments, variables passed into the subroutine
  INTEGER(C_INT), INTENT(IN)          :: i               ! Observation index
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: funvalueI        ! The final computed result
  INTEGER, INTENT(OUT)                :: exitstatus      ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: relerr          ! The final computed result and relative error
  INTEGER, INTENT(IN)                 :: verbose         ! verbosity flag
  INTEGER(C_INT), INTENT(OUT)         :: count_Integration_Regions ! Num int regions

  INTERFACE
    FUNCTION findKmaxSP(j)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE  
      REAL(KIND=C_DOUBLE)   :: findKmaxSP
      INTEGER, INTENT(IN)   :: j
    END FUNCTION findKmaxSP


    SUBROUTINE findKmax(j, kmax, tmax, mmax, mfirst, startTKmax, kmaxL, kmaxR)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER, INTENT(IN)               :: j
      INTEGER(C_INT), INTENT(OUT)       :: mfirst, mmax
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: kmax, tmax
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: startTKmax, kmaxL, kmaxR
    END SUBROUTINE findKmax


    SUBROUTINE improveKmaxSPBounds(i, startx, xL, xR)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: startx
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: xL, xR
    END SUBROUTINE improveKmaxSPBounds


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
      INTEGER(C_INT), INTENT(IN)        :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: xvec(200), wvec(200)
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: Mmatrix(2, 200), Nmatrix(2, 200)
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: West
    END SUBROUTINE acceleratenew


    SUBROUTINE findExactZeros(i, m, tL, tR, zeroSP, zero)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i, m
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: tL, tR, zeroSP
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: zero
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


  END INTERFACE
  ! --- END INTERFACES ---

  ! --- Subroutine Outputs/Inputs ---
  REAL(KIND=C_DOUBLE)             :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200)
  REAL(KIND=C_DOUBLE), VOLATILE   :: West, Wold, Wold2
  REAL(KIND=C_DOUBLE)             :: zeroBoundR, zeroBoundL, zeroStartPoint

  
  ! --- Local Variables ---
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR, zero, aimrerr, epsilon
  INTEGER               :: m, mOld, mmax, mfirst, accMax, minAccRegions
  INTEGER               :: leftOfMax, flip, convergence, stopPreAccelerate
  REAL(KIND=C_DOUBLE)   :: pi
  INTEGER               :: itsacceleration, itsPreAcc

  REAL(KIND=C_DOUBLE)   :: kmax, tmax, startTKmax
  REAL(KIND=C_DOUBLE)   :: df, f, finalTP, front, kmaxL, kmaxR
  REAL(KIND=C_DOUBLE)   :: areaT, area0, area1, areaA, sum, psi
  REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi ! Still using KIND=C_DOUBLE for internal module array access

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
 
  ! --- Initialization ---
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  aimrerr = 1.0E-12_C_DOUBLE
  mOld = 0
  m = 0
  exitstatus = 0
  relerr = 1.0_C_DOUBLE
  convergence = 0
  epsilon = 1.0-13_C_DOUBLE
  mmax = 0

  IF (verbose .EQ. 1) CALL DBLEPR(" FOR 1 < p < 2: p = ", -1, Cp, 1)

  ! FIND kmax, tmax, mmax
  IF (current_y .GE. current_mu) THEN
    ! ************** y > MU   **************
    ! Im k(t) heads down immediately

    IF (verbose .EQ. 1) WRITE(*,*) "** y >= mu"
    kmax = 0.0_C_DOUBLE
    tmax = 0.0_C_DOUBLE
    mmax = 0
    mfirst = -1
    mOld = 0

    IF (verbose .EQ. 1) CALL DBLEPR("** Im k(t) heads down immediately", -1, current_y, 0)
        
    zeroStartPoint = pi / current_y
    leftOfMax = 0
  ELSE
    ! ************** y < MU   **************
    IF (verbose .EQ. 1) CALL DBLEPR("** y < mu: y  = ", -1, current_y, 1)
      
    startTKMax = findKmaxSP(i)
    IF (verbose .EQ. 1) CALL DBLEPR("Find kmax, starting at: ", -1, StartTKmax, 1)
    
    ! For smallp, sometimes very important to have a good starting point and bounds.
    CALL improveKmaxSPBounds(i, startTKmax, kmaxL, kmaxR)
    startTKmax =  (kmaxL + kmaxR) / 2.0_C_DOUBLE
      
  IF (verbose .EQ. 1) CALL DBLEPR("Find kmax, starting (revised) at: ", -1, StartTKmax, 1)
    CALL findKmax(i, kmax, tmax, mmax, mfirst, startTKmax, kmaxL, kmaxR)
  IF (verbose .EQ. 1) CALL DBLEPR("Found kmax: ", -1, kmax, 1)

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
      CALL DBLEPR("** Found(b): kmax =", -1, kmax, 1)
      CALL DBLEPR("             tmax =", -1, tmax, 1)
      CALL INTPR( "             mmax =", -1, mmax, 1)
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

  ! --- Integration initialization ---
  area0 = 0.0_C_DOUBLE
  area1 = 0.0_C_DOUBLE
  areaA = 0.0_C_DOUBLE
  m = mfirst
  ! Find the final turning point of Im/Re k, and start accelerating thereafter      
  WRITE(*,*) "IS this about TPs correct??"
  finalTP = 0.0_C_DOUBLE
  IF ( current_y .LT. current_mu ) THEN
    ! Find where the oscillations settle down, to use acceleration thereafter
    CALL findAccelStart(i, finalTP)
  END IF

  zeroStartPoint = pi / current_y
  ! TRY A NEW ONE!
  front = current_mu ** (1.0_C_DOUBLE - Cp) / ( current_phi * (1.0_C_DOUBLE - Cp))
  zeroStartPoint = front * DTAN( pi * ( 1.0_C_DOUBLE - Cp) / Cp )


  ! 1. INTEGRATE FIRST REGION: area0
  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "*******************************" 
    WRITE(*,*) "1. INTEGRATE: the INITIAL region"
  END IF 
  zeroBoundL = tmax
  zeroBoundR = zeroStartPoint + 3.0_C_DOUBLE * pi / current_y
  IF (verbose .EQ. 1) WRITE(*,*) " Bounds zero; ", zeroBoundL, zeroBoundR

  CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)
  zeroL = 0.0_C_DOUBLE
  zeroR = zero

  ! Find the right-side zero
  CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)

  zeroL =  0.0_C_DOUBLE
  zeroR = zero
  CALL DFgaussq(i, zeroL, zeroR, area0)
  count_Integration_Regions = count_Integration_Regions + 1

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
  area1 = 0.0_C_DOUBLE
  CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

  stopPreAccelerate = 0

  DO WHILE ( stopPreAccelerate .EQ. 0 )
    itsPreAcc = itsPreAcc + 1
    zeroL = zeroR

    zeroStartPoint = (itsPreAcc + 1) * pi / current_y
    zeroBoundL = zeroR
    zeroBoundR = zeroStartPoint + 0.75_C_DOUBLE * pi / current_y

    CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)

    zeroR = zero
    CALL DFgaussq(i, zeroL, zeroR, sum)
    count_Integration_Regions = count_Integration_Regions + 1
    area1 = area1 + sum
        
    IF (verbose .EQ. 1) THEN
      WRITE(*,*) "--- Integrate (m = ", m, ") between " 
      WRITE(*,*) zeroL, "and ", zeroR, "; sum: ", sum
    END IF

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
  West = 3.0_C_DOUBLE
  Wold = 2.0_C_DOUBLE
  Wold2 = 1.0_C_DOUBLE
  itsAcceleration = 0
  areaA = 0.0_C_DOUBLE
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
      zeroR = zeroR * 20.0_C_DOUBLE
    ELSE
      IF (flip .EQ. 1) THEN
        ! FLIPPING to other side of tmax
        zeroStartPoint = tmax + (tmax - zero)
        ! That is, start of the other side of tmax
        zeroL = zero
        zeroR = zeroStartPoint * 20.0_C_DOUBLE
      ELSE
        zeroStartPoint = zeroR
        zeroL = zeroR
        zeroR = zeroR * 10.0_C_DOUBLE
      END IF
    END IF

    CALL findExactZeros(i, m, zeroL, zeroR, zeroStartPoint, zero)
    WRITE(*,*) " What is next line for...???"
    CALL findZeroSmallp(i, zero, f, df)
    
    zeroR = zero
    xvec(itsAcceleration + 1) = zeroR
    IF (verbose .EQ. 1) WRITE(*,*) "  - Integrate between:", zeroL, zeroR
    IF (verbose .EQ. 1) WRITE(*,*) "    : where m = ", m

    CALL DFgaussq(i, zeroL, zeroR, psi)
    count_Integration_Regions = count_Integration_Regions + 1
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
    CALL DBLEPR("* Initial area0: :", -1, area0, 1)
    CALL DBLEPR("* Pre-acc area1: :", -1, area1, 1)
    CALL DBLEPR("*      Acc West: :", -1, West, 1)
    CALL DBLEPR("***       TOTAL: :", -1, areaT, 1)
  END IF
   
  ! WHAT TO DO with relerrr? Might have three rel errors: from initila, pre-acc, acc?
  ! Take largest of the three? ADD?
  ! Assume the argest relative error comes from the acceleration.

  ! So the value returned by the integration  
  funvalueI = -areaT/pi + 0.50_C_DOUBLE
  IF (verbose .EQ. 1) CALL DBLEPR("***    Fun. value:", -1, funvalueI, 1)
    
END SUBROUTINE DFsmallp
