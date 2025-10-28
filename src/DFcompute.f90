SUBROUTINE DFcompute(i, funvalueI, exitstatus, relerr, count_Integration_Regions) 
  USE DFintegrand_MOD, ONLY: DFintegrand
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE rprintf_mod
  
  IMPLICIT NONE
  
 ! --- Dummy Arguments, variables passed into the subroutine
  INTEGER(C_INT), INTENT(IN)        :: i              ! Observation index
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalueI      ! Computed result
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr         ! Estimate of relative error
  INTEGER(C_INT), INTENT(OUT)       :: count_Integration_Regions ! Num int regions
   ! --- INTERFACES: All C-bound routines called by DFbigp:
  INTERFACE
    SUBROUTINE findKmax(j, kmax, tmax, mmax, mfirst, leftOfMax)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER, INTENT(IN)               :: j
      INTEGER(C_INT), INTENT(OUT)       :: mfirst, mmax, leftOfMax
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: kmax, tmax
    END SUBROUTINE findKmax


    SUBROUTINE improveKZeroBounds(i, m, leftOfMax, startx, xL, xR)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i, m, leftOfMax
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: startx
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: xL, xR
    END SUBROUTINE improveKZeroBounds


    SUBROUTINE advanceM(j, m_index, mmax, mOld, leftOfMax, flip) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)             :: j
      INTEGER(C_INT), INTENT(IN)      :: mmax
      INTEGER(C_INT), INTENT(INOUT)   :: m_index
      INTEGER(C_INT), INTENT(OUT)     :: mOld
      INTEGER(C_INT), INTENT(INOUT)   :: leftOfMax
      INTEGER(C_INT), INTENT(OUT)     :: flip
    END SUBROUTINE advanceM
      

    SUBROUTINE DFgaussq(i, a, b, area)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(KIND=C_DOUBLE), INTENT(OUT)      :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)       :: a, b
    END SUBROUTINE DFgaussq


    SUBROUTINE accelerate(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)               :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: xvec(200), wvec(200), Mmatrix(2, 200), Nmatrix(2, 200)
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: West
    END SUBROUTINE accelerate


    SUBROUTINE findExactZeros(i, m, tL, tR, zeroSP, zero, leftOfMax)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)                 :: i, m, leftOfMax
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: tL, tR, zeroSP
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zero
    END SUBROUTINE findExactZeros

  END INTERFACE
  ! --- END INTERFACES ---

  ! Local Variables: All local variables defined here
  INTEGER(C_INT)  :: mmax, mfirst, mOld, accMax
  INTEGER         :: its_Acceleration, count_PreAcc_regions, m, min_Acc_Regions
  INTEGER         :: leftOfMax, flip, convergence_Acc, stop_PreAccelerate
  
  REAL(KIND=C_DOUBLE)   :: kmax, tmax, aimrerr
  REAL(KIND=C_DOUBLE)   :: epsilon, areaT, pi, psi, zero, t_Start_Point
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR, area0, area1, sumA
  REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE)   :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200)
  REAL(KIND=C_DOUBLE)   :: West, Wold, Wold2
  REAL(KIND=C_DOUBLE)   :: zeroBoundR, zeroBoundL, zeroStartPoint
  

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
  convergence_Acc = 0
  epsilon = 1.0E-12_C_DOUBLE
  Mmatrix = 0.0_C_DOUBLE
  Nmatrix = 0.0_C_DOUBLE
  xvec = 0.0_C_DOUBLE
  wvec = 0.0_C_DOUBLE
  mmax = 0
  zeroStartPoint = 0.0_C_DOUBLE
  
  IF (Cverbose) THEN
    CALL DBLEPR("********* Computing for p =", -1, Cp, 1)
  END IF
  
  ! --- Find kmax, tmax, mmax ---
  CALL findKmax(i, kmax, tmax, mmax, mfirst, leftOfMax)

  IF (Cverbose) THEN
    CALL DBLEPR("  -        kmax:", -1, kmax, 1 )
    CALL DBLEPR("  -        tmax:", -1, tmax, 1 )
    CALL INTPR( "  -        mmax:", -1, mmax, 1 )
    CALL INTPR( "  - left of max:", -1, leftOfMax, 1 )
  END IF
  

  ! INTEGRATION
  ! There are three integration regions:
  !   1. The *initial* area, which is not between zeros of Im{k(t)}: area0
  !   2. The initial area *before* Sidi acceleration is invoked: area1
  !      (For instance, wait until Im{k(t)} is on the downturn.)
  !   3. The area thereafter, upon which Sidi acceleration is
  !      applied; the area returned by acceleration is West

  ! --- Integration initialization ---
  area0 = 0.0_C_DOUBLE
  area1 = 0.0_C_DOUBLE
  zero  = 0.0_C_DOUBLE
  

  ! ---------------------------------------------------------
  ! --- 1. INTEGRATE FIRST REGION: area0 ---
  count_Integration_Regions = 1

  ! Find starting point for the first zero
  m = mfirst
  t_Start_Point = tmax + pi / current_y  
  IF (leftOfMax .EQ. 1) THEN
    zeroBoundL = 0.0_C_DOUBLE
    zeroBoundR = t_Start_Point * 2.0_C_DOUBLE
  ELSE
    ! Searching to the right of tmax
    zeroBoundL = tmax
    zeroBoundR = t_Start_Point * 2.0_C_DOUBLE
  END IF

  ! Find the zero
  zeroL = 0.0_C_DOUBLE
  CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, t_Start_Point, zeroR, leftOfMax)

  ! Find the area
  CALL DFgaussq(i, zeroL, zeroR, area0)

  IF (Cverbose) THEN
    CALL DBLEPR("  *** INITIAL area:", -1, area0, 1 )
    CALL DBLEPR("        between t =", -1, zeroL, 1 )
    CALL DBLEPR("            and t =", -1, zeroR, 1 )
    CALL INTPR( "          using m =", -1, m, 1)
  END IF


  ! ---------------------------------------------------------
  ! --- 2. INTEGRATE: the PRE-ACCELERATION regions: area1 ---
  count_PreAcc_regions = 0
  IF (mfirst .EQ. -1) THEN
    ! count_PreAcc_regions = count_PreAcc_regions + 1
    area1 = 0.0_C_DOUBLE
    mOld = m

    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
    IF (Cverbose) CALL DBLEPR("  - No PRE-ACC area for y:", -1, current_y, 1 )

  ELSE
    area1 = 0.0_C_DOUBLE
    mOld = m

    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

    stop_PreAccelerate = 0
    DO WHILE (stop_PreAccelerate .EQ. 0)
      count_PreAcc_regions = count_PreAcc_regions + 1

      IF (leftOfMax .EQ. 1) THEN
        zeroBoundL = zeroR
        zeroBoundR = zeroR * 10.0_C_DOUBLE
      ELSE
        zeroBoundL = tmax
        zeroBoundR = zeroR * 20.0_C_DOUBLE
      END IF

      zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
      zeroL = zeroR

!ONLY NEEDED FOR SMALL P
!CALL improveKZeroBounds(i, m, leftOfMax, zeroStartPoint, zeroBoundL, zeroBoundR)
      CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero, leftOfMax)

      zeroR = zero

      CALL DFgaussq(i, zeroL, zeroR, sumA)
      area1 = area1 + sumA
      count_Integration_Regions = count_Integration_Regions + 1
      IF (Cverbose) THEN
        CALL DBLEPR("  *** PRE-ACC area:", -1, sumA, 1 )
        CALL DBLEPR("        between t =", -1, zeroL, 1 )
        CALL DBLEPR("            and t =", -1, zeroR, 1 )
        CALL INTPR( "          using m =", -1, m, 1 )
      END IF

      ! Stop condition for pre-acceleration.
      ! Ensure that we have passed the peak of Im k(t), so that acceleration can be used
      IF ( (count_PreAcc_regions .GE. 2) .AND. (zeroL .GT. tmax) ) stop_PreAccelerate = 1
      
      mOld = m
      CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
    
    END DO 
  END IF


  ! ----------------------------------------------------
  ! --- 3. INTEGRATE: the ACCELERATION regions: West ---
  IF (Cverbose) CALL DBLEPR(" - ACCELERATING for y:", -1, current_y, 1)

  ! Initialisation of acceleration 
  West = 3.0_C_DOUBLE
  Wold = 2.0_C_DOUBLE
  Wold2 = 1.0_C_DOUBLE
  its_Acceleration = 0
  convergence_Acc = 0
  accMax = 40           ! Maximum number of regions in acceleration; arbitrary
  min_Acc_Regions = 3     ! Minimum number of acceleration regions to use
  
  ! This will be the very first, left-most value of t used in acceleration
  xvec(1) = zeroR 

  DO WHILE (convergence_Acc .EQ. 0)
    its_Acceleration = its_Acceleration + 1

    zeroL = zeroR

    IF (flip .EQ. 1) THEN
      ! FLIPPING to other side of tmax
      IF (Cverbose) CALL DBLEPR("  - Flipping to other side of tmax", -1, tmax, 1)
      zeroStartPoint = tmax + (tmax - zeroR)
      ! That is, start of the other side of tmax
      zeroBoundL = zeroR 
      zeroBoundR = zeroStartPoint * 2.0_C_DOUBLE
    ELSE
      zeroStartPoint = zeroR 
      zeroBoundL = zeroR * 0.99_C_DOUBLE
      zeroBoundR = zeroR * 2.0_C_DOUBLE
    END IF

    ! Find the exact zero
    CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero, leftOfMax)

    zeroR = zero

    xvec(its_Acceleration + 1) = zeroR
    
    IF (Cverbose) CALL INTPR("  - m =:", -1, m, 1 )
    
    CALL DFgaussq(i, zeroL, zeroR, psi)
    count_Integration_Regions = count_Integration_Regions + 1
    ! psi: area of the latest region

    wvec(its_Acceleration) = psi
    Wold2 = Wold
    Wold = West
    CALL accelerate(xvec, wvec, its_Acceleration, Mmatrix, Nmatrix, West)

    ! Check for convergence_Acc
    relerr = (DABS(West - Wold) + DABS(West - Wold2)) / (DABS(West) + epsilon)
    IF (Cverbose) THEN
      CALL DBLEPR("  *** Acceleration area:", -1, psi, 1)
      CALL DBLEPR("             between t =", -1, zeroL, 1)
      CALL DBLEPR("                 and t =", -1, zeroR, 1)
      CALL INTPR( "               using m =", -1, m, 1)
      CALL DBLEPR("          with rel err =", -1, relerr, 1)
    END IF

    ! Declare convergence_Acc of we have sufficient regions, and relerr estimate is small
    IF ( (its_Acceleration .GE. min_Acc_Regions) .AND. &
         (relerr .LT. aimrerr) ) THEN
      IF (Cverbose) CALL DBLEPR(   "         rel err =:", -1, relerr, 1)
      convergence_Acc = 1
    END IF
    
!    IF (its_Acceleration .GT. 2) STOP

    mOld = m
    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

    ! NOTE: If convergence_Acc is NOT TRUE, the loop continues.
  END DO  

  areaT = area0 + area1 + West
  
  IF (Cverbose) THEN
    CALL DBLEPR("* Initial area0: :", -1, area0, 1)
    CALL DBLEPR("* Pre-acc area1: :", -1, area1, 1)
    CALL DBLEPR("*      Acc West: :", -1, West, 1)
    CALL DBLEPR("***       TOTAL: :", -1, areaT, 1)
  END IF

  ! We have the value of the integral in the PDF/CDF calculation.
  ! So now work out the actualy PDF/CDF
  
  IF (Cpdf) THEN
    funvalueI = areaT/pi 
  ELSE
    funvalueI = -areaT/pi + 0.5E0_C_DOUBLE
  END IF  
  IF (Cverbose) CALL DBLEPR("***    Fun. value:", -1, funvalueI, 1)

  RETURN

END SUBROUTINE DFcompute

