SUBROUTINE DFbigp(i, funvalueI, exitstatus, relerr, verbose, count_Integration_Regions) 
  USE DFintegrand_MOD, ONLY: DFintegrand
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE rprintf_mod
  
  IMPLICIT NONE
  
 ! --- Dummy Arguments, variables passed into the subroutine
  INTEGER(C_INT), INTENT(IN)        :: i              ! Observation index
  INTEGER(C_INT), INTENT(IN)        :: verbose        ! Assuming IN for verbosity flag
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalueI      ! Computed result
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr         ! Estimate of relative error
  INTEGER(C_INT), INTENT(OUT)       :: count_Integration_Regions ! Num int regions
   ! --- INTERFACES: All C-bound routines called by DFbigp:
  INTERFACE
    FUNCTION findKmaxSP(j) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE  
      REAL(KIND=C_DOUBLE)          :: findKmaxSP
      INTEGER, INTENT(IN)   :: j
    END FUNCTION findKmaxSP


    SUBROUTINE findKmax(j, kmax, tmax, mmax, mfirst, startTKmax) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER, INTENT(IN)               :: j
      INTEGER(C_INT), INTENT(OUT)       :: mfirst, mmax
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: kmax, tmax
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: startTKmax
    END SUBROUTINE findKmax

      ! 3. Subroutine to advance the iteration index m
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
      REAL(KIND=C_DOUBLE)             :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)        :: a, b
    END SUBROUTINE DFgaussq


    SUBROUTINE acceleratenew(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)               :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: xvec(200), wvec(200), Mmatrix(2, 200), Nmatrix(2, 200)
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: West
    END SUBROUTINE acceleratenew


    SUBROUTINE findExactZeros(i, m, tL, tR, zeroL, zeroR)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)         :: i, m
      REAL(KIND=C_DOUBLE), INTENT(IN)    :: tL, tR
      REAL(KIND=C_DOUBLE), INTENT(OUT)   :: zeroL, zeroR
    END SUBROUTINE findExactZeros

  END INTERFACE
  ! --- END INTERFACES ---

  ! Local Variables: All local variables defined here
  INTEGER(C_INT)  :: mmax, mfirst, mOld, accMax
  INTEGER         :: itsAcceleration, itsPreAcc, m, minAccRegions
  INTEGER         :: leftOfMax, flip, convergence, stopPreAccelerate
  
  REAL(KIND=C_DOUBLE) :: kmax, startTKmax, tmax, aimrerr
  REAL(KIND=C_DOUBLE) :: epsilon, areaT, pi, psi, zero
  REAL(KIND=C_DOUBLE) :: zeroL, zeroR, area0, area1, sumA
  REAL(KIND=C_DOUBLE) :: current_y, current_mu, current_phi ! Can still using KIND=C_DOUBLE for internal module array access
  REAL(KIND=C_DOUBLE)              :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200)
  REAL(KIND=C_DOUBLE)    :: West, Wold, Wold2
  REAL(KIND=C_DOUBLE)    :: zeroBoundR, zeroBoundL, zeroStartPoint
  

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
  epsilon = 1.0E-12_C_DOUBLE
  Mmatrix = 0.0_C_DOUBLE
  Nmatrix = 0.0_C_DOUBLE
  xvec = 0.0_C_DOUBLE
  wvec = 0.0_C_DOUBLE
  mmax =- 0
  zeroStartPoint = 0.0_C_DOUBLE
 
  IF (verbose .EQ. 1) THEN
    CALL DBLEPR("********* Computing for p =", -1, Cp, 1)
  END IF
  
  ! --- Find kmax, tmax, mmax ---
!write(*,*) "DEBUG: Index i before advanceM is: ", i
  IF (current_y .GE. current_mu) THEN
    IF (verbose .EQ. 1) THEN
      CALL DBLEPR("  - with y >= mu: y =", -1, current_y, 1 )
    END IF
    
    kmax = 0.0_C_DOUBLE
    tmax = 0.0_C_DOUBLE
    mmax = 0
    mfirst = -1
    mOld = 0

    zeroStartPoint = pi / current_y
    leftOfMax = 0
  ELSE
    IF (verbose .EQ. 1) THEN
      CALL DBLEPR("  - with < mu: y =", -1, current_y, 1 )
    END IF
    
    startTKmax = findKmaxSP(i)
!  WRITE(*,*) "!!! TRACE 1: zeroStartPoint AFTER findKmax:", zeroStartPoint
 
    CALL findKmax(i, kmax, tmax, mmax, mfirst, startTKmax)

    IF (verbose .EQ. 1) THEN
      CALL DBLEPR("  - kmax:", -1, kmax, 1 )
      CALL DBLEPR("  - tmax:", -1, tmax, 1 )
      CALL INTPR( "  - mmax:", -1, mmax, 1 )
    END IF

    leftOfMax = 1
    IF (mmax .EQ. 0) THEN
      mfirst = 0
      mOld = 0
      zeroStartPoint = tmax + pi/current_y
      zeroStartPoint = zeroStartPoint
      leftOfMax = 0
    ELSE
      mfirst = 1
      mOld = 0
      zeroStartPoint = pi / (current_mu - current_y)
      zeroStartPoint = zeroStartPoint
      mOld = m
  
      ! CRITICAL: advanceM must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
      CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
  
    END IF
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

  ! --- 1. INTEGRATE FIRST REGION: area0 ---
  count_Integration_Regions = 1
  zeroBoundL = tmax ! 0.0_C_DOUBLE  ! Should be tmax???
  zeroBoundR = zeroStartPoint * 2.0_C_DOUBLE
  m = mfirst

  CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)

  zeroL = 0.0_C_DOUBLE
  zeroR = zero

  CALL DFgaussq(i, zeroL, zeroR, area0)
  IF (verbose .EQ. 1) THEN
    CALL DBLEPR("  * INITIAL area:", -1, area0, 1 )
    CALL DBLEPR("      between t =", -1, zeroL, 1 )
    CALL DBLEPR("          and t =", -1, zeroR, 1 )
    CALL INTPR( "        using m =", -1, m, 1)
  END IF

  ! --- 2. INTEGRATE: the PRE-ACCELERATION regions: area1 ---
  itsPreAcc = 0
  IF (mfirst .EQ. -1) THEN
    itsPreAcc = itsPreAcc + 1
    area1 = 0.0_C_DOUBLE
    mOld = m

    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
    IF (verbose .EQ. 1) CALL DBLEPR("  - No PRE-ACC area for y:", -1, current_y, 1 )

  ELSE
    area1 = 0.0_C_DOUBLE
    mOld = m

    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

    stopPreAccelerate = 0
    DO WHILE (stopPreAccelerate .EQ. 0)
      itsPreAcc = itsPreAcc + 1

      IF (leftOfMax .EQ. 1) THEN
        zeroBoundL = zeroR
        zeroBoundR = zeroR * 10.0_C_DOUBLE
      ELSE
        zeroBoundL = tmax
        zeroBoundR = zeroR * 20.0_C_DOUBLE
      END IF
      zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
      
      zeroL = zeroR
  
      CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)

      zeroR = zero

      CALL DFgaussq(i, zeroL, zeroR, sumA)
      area1 = area1 + sumA
      count_Integration_Regions = count_Integration_Regions + 1
      IF (verbose .EQ. 1) THEN
        CALL DBLEPR("  - PRE-ACC area:", -1, sumA, 1 )
        CALL DBLEPR("      between t =", -1, zeroL, 1 )
        CALL DBLEPR("          and t =", -1, zeroR, 1 )
        CALL INTPR( "        using m =", -1, m, 1 )
      END IF

      ! Stop condition for pre-acceleration.
      IF (itsPreAcc .GE. 2) stopPreAccelerate = 1
      
      mOld = m
      CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
    
    END DO 
  END IF

  ! --- 3. INTEGRATE: the ACCELERATION regions: West ---
  IF (verbose .EQ. 1) CALL rprintf_double(" - ACCELERATING for y:", current_y)

  ! Initialisation
  West = 3.0_C_DOUBLE
  Wold = 2.0_C_DOUBLE
  Wold2 = 1.0_C_DOUBLE
  itsAcceleration = 0
  convergence = 0
  accMax = 40           ! Maximum number of regions in accelerationl arbitrary
  minAccRegions = 3     ! Minimum number of acceleration regions to use
  
  ! This will be the very first, left-most value of t used in acceleration
  xvec(1) = zeroR 

  DO WHILE (convergence .EQ. 0)

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

    ! Find the first exact zero
    CALL findExactZeros(i, m, zeroL, zeroR, zeroStartPoint, zero)

    IF (leftOfMax .EQ. 1) THEN
      zeroR = zero
    ELSE
      zeroR = zero
    END IF

    xvec(itsAcceleration + 1) = zeroR
    
    IF (verbose .EQ. 1) THEN
      CALL INTPR("  - m =:", -1, m, 1 )
    END IF

    CALL DFgaussq(i, zeroL, zeroR, psi)
    count_Integration_Regions = count_Integration_Regions + 1
    ! psi: area of the latest region
    wvec(itsAcceleration) = psi

    Wold2 = Wold
    Wold = West
    CALL acceleratenew(xvec, wvec, itsAcceleration, Mmatrix, Nmatrix, West)

    ! Check for convergence
    relerr = (DABS(West - Wold) + DABS(West - Wold2)) / (DABS(West) + epsilon)
    IF (verbose .EQ. 1) THEN
      CALL DBLEPR("  - Acceleration area:", -1, psi, 1)
      CALL DBLEPR("           between t =", -1, zeroL, 1)
      CALL DBLEPR("               and t =", -1, zeroR, 1)
      CALL INTPR( "             using m =", -1, m, 1)
      CALL DBLEPR("        with rel err =", -1, relerr, 1)
    END IF

    ! Declare convergence of we have sufficient regions, and relerr estimate is small
    IF ( (itsAcceleration .GE. minAccRegions) .AND. &
         (relerr .LT. aimrerr) ) THEN
      IF (verbose .EQ. 1) CALL rprintf_double(   "         rel err =:", relerr)
      convergence = 1
    END IF
  
    mOld = m
    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
    
    ! NOTE: If convergence is NOT TRUE, the loop continues.
  END DO  

  areaT = area0 + area1 + West
  
  IF (verbose .EQ. 1) THEN
    CALL DBLEPR("* Initial area0: :", -1, area0, 1)
    CALL DBLEPR("* Pre-acc area1: :", -1, area1, 1)
    CALL DBLEPR("*      Acc West: :", -1, West, 1)
    CALL DBLEPR("***       TOTAL: :", -1, areaT, 1)
  END IF

  ! We have the value of the integral in the CDF calculation.
  ! So now work out the CDF
  funvalueI = (-1.0_C_DOUBLE/pi) * areaT + 0.5_C_DOUBLE
    
  IF (verbose .EQ. 1) CALL DBLEPR("***    Fun. value:", -1, funvalueI, 1)

  RETURN

END SUBROUTINE DFbigp

