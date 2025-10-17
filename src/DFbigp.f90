SUBROUTINE DFbigp(i, funvalue, exitstatus, relerr, verbose) 
  USE DFintegrand_MOD, ONLY: DFintegrand
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
  IMPLICIT NONE
  
 ! --- Dummy Arguments, variables passed into the subroutine
  INTEGER(C_INT), INTENT(IN)                      :: i              ! Observation index
  INTEGER(C_INT), INTENT(IN)                      :: verbose        ! Assuming IN for verbosity flag
  INTEGER(C_INT), INTENT(OUT)                     :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), DIMENSION(*), INTENT(OUT)  :: funvalue       ! Computed result
  REAL(KIND=C_DOUBLE), INTENT(OUT)                :: relerr         ! Estimate of relative error

   ! --- INTERFACES: All C-bound routines called by DFbigp:
  INTERFACE
    FUNCTION findKmaxSP(j) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE  
      REAL(KIND=8)          :: findKmaxSP
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
      REAL(KIND=8), INTENT(IN)        :: a, b
    END SUBROUTINE DFgaussq


    SUBROUTINE acceleratenew(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      INTEGER, INTENT(IN)      :: nzeros
      REAL(KIND=8), INTENT(IN) :: xvec(200), wvec(200), Mmatrix(2, 200), Nmatrix(2, 200), West
    END SUBROUTINE acceleratenew


    SUBROUTINE findExactZeros(i, m, tL, tR, zeroL, zeroR)
      INTEGER, INTENT(IN)         :: i, m
      REAL(KIND=8), INTENT(IN)    :: tL, tR
      REAL(KIND=8), INTENT(OUT)   :: zeroL, zeroR
    END SUBROUTINE findExactZeros

  END INTERFACE
  ! --- END INTERFACES ---

  ! Local Variables: All local variables defined here
  INTEGER(C_INT)  :: mmax, mfirst, mOld, accMax
  INTEGER         :: itsAcceleration, itsPreAcc, m, minAccRegions
  INTEGER         :: leftOfMax, flip, convergence, stopPreAccelerate
  
  REAL(KIND=8) :: kmax, startTKmax, tmax, aimrerr
  REAL(KIND=8) :: epsilon, areaT, pi, psi, zero
  REAL(KIND=8) :: zeroL, zeroR, area0, area1, areaA, sum
  REAL(KIND=8) :: current_y, current_mu, current_phi ! Can still using KIND=8 for internal module array access
  REAL(KIND=8)              :: Mmatrix(2, 200), Nmatrix(2, 200), xvec(200), wvec(200)
  REAL(KIND=8), VOLATILE    :: West, Wold, Wold2
  REAL(KIND=8)              :: zeroBoundR, zeroBoundL, zeroStartPoint
  

  ! --- Initialization ---
  pi = 4.0D0 * DATAN(1.0D0)
  aimrerr = 1.0D-12
  mOld = 0
  m = 0

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
 
  IF (verbose .EQ. 1) WRITE(*,*) " FOR p > 2"
    
    exitstatus = 0
    relerr = 1.0_8
    convergence = 0
    epsilon = 1.0d-12
    
    ! --- Find kmax, tmax, mmax ---
    IF (Cy(i) .GE. Cmu(i)) THEN
      IF (verbose .EQ. 1) WRITE(*,*) "** y >= mu"
      
      kmax = 0.0_8
      tmax = 0.0_8
      mmax = 0
      mfirst = -1
      mOld = 0
      IF (verbose .EQ. 1) WRITE(*,*) "** Im k(t) heads down immediately"
      
      zeroStartPoint = pi / current_y
      leftOfMax = 0
    ELSE
      IF (verbose .EQ. 1) WRITE(*,*) "** y < mu"
    
      ! CRITICAL: findKmaxSP must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
      startTKmax = findKmaxSP(i)

      IF (verbose .EQ. 1) WRITE(*,*) "Starting t for finding kmax: ", startTKmax
      CALL findKmax(i, kmax, tmax, mmax, mfirst, startTKmax)
    
      IF (verbose .EQ. 1) THEN
        WRITE(*,*) "** Found(b): kmax =", kmax
        WRITE(*,*) "             tmax =", tmax
        WRITE(*,*) "             mmax =", mmax
      END IF
    
      leftOfMax = 1
      IF (mmax .EQ. 0) THEN
        mfirst = 0
        mOld = 0
        zeroStartPoint = tmax + pi/Cy(i)
        leftOfMax = 0
      ELSE
        mfirst = 1
        mOld = 0
        zeroStartPoint = pi / (current_mu -current_y)
        mOld = m
  
        ! CRITICAL: advanceM must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
        CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
  
      END IF
  END IF
  
  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "            mfirst =", mfirst
    WRITE(*,*) "              StPt =", zeroStartPoint
    WRITE(*,*) "--- (Deal with returned errors, non-convergence)"
  END IF
  
  ! --- Integration initialization ---
  area0 = 0.0_8
  area1 = 0.0_8
  areaA = 0.0_8

  ! --- 1. INTEGRATE FIRST REGION: area0 ---
  IF (verbose .EQ. 1) WRITE(*,*) "*******************************"
  IF (verbose .EQ. 1) WRITE(*,*) "1. INTEGRATE: the INITIAL region"
  
  zeroBoundL = 0.0_8
  zeroBoundR = zeroStartPoint * 2.0_8

  m = mfirst ! This line caused the error; now fixed by INOUT
  ! CRITICAL: findExactZeros must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
!  WRITE(*,*) "DFbigp: ABOUT TO CALL findExactZeros"
  CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)
!  WRITE(*,*) "DFbigp: END CALL findExactZeros"
  zeroL = 0.0_8
  zeroR = zero

  IF (verbose .EQ. 1) WRITE(*,*) "  - Between ", zeroL, zeroR
  ! CRITICAL: DFintegrand must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
  ! CRITICAL: gaussq must be updated to pass these parameters to DFintegrand
  CALL DFgaussq(i, zeroL, zeroR, area0)
  IF (verbose .EQ. 1) WRITE(*,*) "  - Initial area is", area0
  
  ! --- 2. INTEGRATE: the PRE-ACCELERATION regions: area1 ---
  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "*******************************"
    WRITE(*,*) "2. INTEGRATE: the PRE-ACCELERATION regions"
  END IF
  
  itsPreAcc = 0
  IF (mfirst .EQ. -1) THEN
    itsPreAcc = itsPreAcc + 1
    IF (verbose .EQ. 1) WRITE(*,*) "  > Not using pre-acceleration area"
    area1 = 0.0_8
    mOld = m
    ! CRITICAL: advanceM must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

  ELSE
    area1 = 0.0_8
    mOld = m
    ! CRITICAL: advanceM must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)

    stopPreAccelerate = 0
    
    DO WHILE (stopPreAccelerate .EQ. 0)
      itsPreAcc = itsPreAcc + 1

      IF (leftOfMax .EQ. 1) THEN
        zeroBoundL = zeroR
        zeroBoundR = zeroR * 10.0_8
      ELSE
        zeroBoundL = tmax
        zeroBoundR = zeroR * 20.0_8
      END IF
      zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_8
      
      zeroL = zeroR
      ! CRITICAL: findExactZeros must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
      CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero)

      zeroR = zero
      IF (verbose .EQ. 1) WRITE(*,*) "--- Integrate (m = ", m, ") between "
      
      ! CRITICAL: gaussq must be updated to pass these parameters
      CALL DFgaussq(i, zeroL, zeroR, sum)
      area1 = area1 + sum
      IF (verbose .EQ. 1) WRITE(*,*) zeroL, "and ", zeroR, ": ", sum

      ! Stop condition for pre-acceleration.
      IF (itsPreAcc .GE. 2) stopPreAccelerate = 1
      
      mOld = m
      ! CRITICAL: advanceM must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
      CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
    
    END DO 
  END IF

  ! --- 3. INTEGRATE: the ACCELERATION regions: areaA ---
  IF (verbose .EQ. 1) WRITE(*,*) "*******************************"
  IF (verbose .EQ. 1) WRITE(*,*) "3. INTEGRATE: the ACCELERATION"
  
  ! Initialisation
  West = 3.0_8
  Wold = 2.0_8
  Wold2 = 1.0_8
  itsAcceleration = 0
  areaA = 0.0_8
  convergence = 0
  accMax = 40           ! Maximum number of regions in accelerationl arbitrary
  minAccRegions = 3     ! Minimum number of acceleration regions to use
  
  ! This will be the very first, left-most value of t used in acceleration
  xvec(1) = zeroR 

  DO WHILE (convergence .EQ. 0)
    IF (verbose .EQ. 1) WRITE(*,*) "  --------------------- "
    IF (verbose .EQ. 1) WRITE(*,*) "  Next tail region"

    itsAcceleration = itsAcceleration + 1
    
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

    ! Find the first exact zero
    CALL findExactZeros(i, m, zeroL, zeroR, zeroStartPoint, zero)

    IF (leftOfMax .EQ. 1) THEN
      zeroR = zero
    ELSE
      zeroR = zero
    END IF

    xvec(itsAcceleration + 1) = zeroR
    
    IF (verbose .EQ. 1) WRITE(*,*) "  - Integrate (m = ", m, "):", zeroL, zeroR

    ! CRITICAL: gaussq must be updated to pass these parameters
    CALL DFgaussq(i, zeroL, zeroR, psi)
    ! psi: area of the latest region
    wvec(itsAcceleration) = psi
    IF (verbose .EQ. 1) WRITE(*,*) "  - Area between zeros is:", psi

    WRITE(*,*) "PRE-acc 1: ", West, Wold, Wold2
    Wold2 = Wold
    WRITE(*,*) "PRE-acc 2: ", West, Wold, Wold2
    Wold = West
    WRITE(*,*) "PRE-acc 3: ", West, Wold, Wold2
    CALL acceleratenew(xvec, wvec, itsAcceleration, Mmatrix, Nmatrix, West)
    WRITE(*,*) "POST-acc: updated to: ", West, Wold, Wold2
    
    ! Check for convergence
    relerr = (DABS(West - Wold) + DABS(West - Wold2)) / (DABS(West) + epsilon)

    WRITE(*,*) "Rel error:", relerr
    
    IF (verbose .EQ. 1) THEN 
      WRITE(*,*) "  - iteration", itsAcceleration, ":", West
      WRITE(*,*) "  - Estimate of tail area:", West
    END IF
    
    ! Declare convergemce of we have sufficient regions, and relerr estimate is small
    IF ( (itsAcceleration .GE. minAccRegions) .AND. &
         (relerr .LT. aimrerr) ) THEN
      IF (verbose .EQ. 1) WRITE(*,*) "  Relerr is", relerr
      convergence = 1
    END IF
    WRITE(*,*) "CONVERGENCE", convergence
    mOld = m
    ! CRITICAL: advanceM must accept parameters Cp, Cy, Cmu, Cphi, pSmall, m
    CALL advanceM(i, m, mmax, mOld, leftOfMax, flip)
    
    if (verbose .EQ. 1) WRITE(*,*) "--------------------------------"
    ! NOTE: If convergence is NOT TRUE, the loop continues.
  END DO  

  IF (convergence .EQ. 0) THEN
    IF (verbose .EQ. 1) WRITE(*,*) "  - Computing for p > 2:", Cp
  END IF

  areaA = West
  areaT = area0 + area1 + areaA
  
  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "SUMMARY:"
    WRITE(*,*) "  Area0 ", area0
    WRITE(*,*) "  Area1 ", area1, "(", itsPreAcc, "regions)"
    WRITE(*,*) "  AreaA ", areaA, "(", itsAcceleration, " its)"
    WRITE(*,*) "  TOTAL ", areaT
    WRITE(*,*) "FIX rel err: |A|.relA + ... + |C|.relC/|A+B+C|"
  END IF
  
  ! We have the value of the integral in the CDF calculation.
  ! So now work out the CDF
  funvalue(i) = (-1.0_8/pi) * areaT + 0.5_8
    
  IF (verbose .EQ. 1) THEN
    WRITE(*,*) "FINAL AREA: The cdf value is", funvalue(i)
    WRITE(*,*) "DFbigp: funvalue, exitstatus, relerr"
    WRITE(*,*) funvalue(i), exitstatus, relerr
  END IF

  
  RETURN

END SUBROUTINE DFbigp

