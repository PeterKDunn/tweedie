SUBROUTINE ComputeTwIntegral(i, funvalueI, exitstatus, relerr, count_Integration_Regions) 
  ! Compute the value of the integrals in the Fourier-inversion expressions for the PDF and CDF

  USE Integrands_MOD, ONLY: Integrands
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE rprintf_mod
  
  IMPLICIT NONE
  
  INTEGER(C_INT), INTENT(IN)        :: i              ! Observation index
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalueI      ! Computed result
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr         ! Estimate of relative error
  INTEGER(C_INT), INTENT(OUT)       :: count_Integration_Regions ! Num int regions


  ! Local Variables: All local variables defined here
  INTEGER(C_INT)    :: mmax, mfirst, mOld, accMax
  INTEGER(C_INT)    :: count_PreAcc_regions, count_Acc_Regions
  INTEGER(C_INT)    :: m, min_Acc_Regions
  INTEGER(C_INT)    :: leftOfMax, flip_To_Other_Side, convergence_Acc
  LOGICAL(C_BOOL)   :: converged_Accelerating, converged_Pre
  REAL(KIND=C_DOUBLE)   :: kmax, tmax, aimrerr, tStartAcc
  REAL(KIND=C_DOUBLE)   :: epsilon, areaT, pi, psi, zero
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR, area0, area1, areaA, sumA
  REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE)   :: Mmatrix(2, 501), Nmatrix(2, 501), xvec(501), wvec(501)
  REAL(KIND=C_DOUBLE)   :: TMP
  REAL(KIND=C_DOUBLE)   :: zeroStartPoint
  
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
    

    SUBROUTINE GaussQuadrature(i, a, b, area)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(KIND=C_DOUBLE), INTENT(OUT)      :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)       :: a, b
    END SUBROUTINE GaussQuadrature


    SUBROUTINE accelerate(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)                 :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: xvec(501), wvec(501), Mmatrix(2, 501), Nmatrix(2, 501)
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: West
    END SUBROUTINE accelerate


    SUBROUTINE findExactZeros(i, m, tL, tR, zeroSP, zero, leftOfMax)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)                 :: i, m, leftOfMax
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: tL, tR, zeroSP
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zero
    END SUBROUTINE findExactZeros

  END INTERFACE


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  IF (Cverbose) THEN
    CALL DBLEPR("*** Computing for p =", -1, Cp, 1)
    CALL DBLEPR("*** Computing for y =", -1, current_y, 1)
  END IF


  ! --- Initialization ---
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  aimrerr = 1.0E-10_C_DOUBLE ! 1.0E-12_C_DOUBLE
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
  count_Integration_Regions = 0_C_INT
  mmax = 0
  zeroStartPoint = 0.0_C_DOUBLE

  ! --- Integration initialization ---
  area0 = 0.0_C_DOUBLE
  area1 = 0.0_C_DOUBLE
  zero  = 0.0_C_DOUBLE

  ! --- Find kmax, tmax, mmax ---
  CALL findKmax(i, kmax, tmax, mmax, mfirst, leftOfMax)
  IF (Cverbose) THEN
    CALL DBLEPR("  -            kmax:", -1, kmax, 1 )
    CALL DBLEPR("  -            tmax:", -1, tmax, 1 )
    CALL INTPR( "  -            mmax:", -1, mmax, 1 )
    CALL INTPR( "  -     left of max:", -1, leftOfMax, 1 )
    CALL INTPR( "  - first zero at m:", -1, mfirst, 1 )
  END IF
  
  ! INTEGRATION
  ! Three integration regions:
  !   1. The *initial* area, which is not between zeros of Im{k(t)}: area0
  !   2. The initial area *before* Sidi acceleration is invoked: area1
  !      (For instance, wait until Im{k(t)} is on the downturn.)
  !   3. The area thereafter, upon which Sidi acceleration is
  !      applied; the area returned by acceleration is areaA


  ! ----------------------------------------------------------------------------
  ! --- 1. INTEGRATE FIRST (sometimes non-standard) REGION: area0 ---
  CALL integrateFirstRegion(mfirst, leftOfMax,          & ! INPUTS
                            area0, zeroR)                 ! OUTPUTS
  count_Integration_Regions = 1

  IF (Cverbose) THEN
    CALL DBLEPR("Initial region area:", -1, area0, 1)
    CALL DBLEPR("      between 0 and:", -1, zeroR, 1)
    CALL INTPR( "      using right m:", -1, mfirst, 1)
  END IF



  ! ----------------------------------------------------------------------------
  ! --- 2. INTEGRATE: the PRE-ACCELERATION regions: area1 ---

  ! Find where to start accelerating
  IF (Cp .GT. 2.0_C_DOUBLE) THEN
    tStartAcc = 0.0_C_DOUBLE
  ELSE
    CALL findAccelStart(i, tmax, tStartAcc)
  END IF
  IF (Cverbose) CALL DBLEPR(" Start accelerating after:", -1, tStartAcc, 1)


  zeroL = zeroR  ! The last region's right-side zero is next region's left-side zero
  CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)
  CALL integratePreAccRegions(m, mfirst, leftOfMax, zeroL,  tmax, tStartAcc,        & ! INPUTS
                              area1, zeroR, count_PreAcc_regions,  converged_Pre)     ! OUTPUTS
  count_Integration_Regions = count_Integration_Regions + count_PreAcc_regions

  IF (Cverbose) THEN
    CALL DBLEPR(" Pre-acc area:", -1, area1, 1)
    CALL DBLEPR("      between:", -1, zeroL, 1)
    CALL DBLEPR("          and:", -1, zeroR, 1)
    CALL INTPR( "using right m:", -1, m,     1)
  
    IF (converged_Pre) THEN
      CALL DBLEPR(" Accelerating not needed; convergence by t =", -1, zeroR, 1)
    ELSE
      CALL DBLEPR(" Accelerating starting after t =", -1, zeroR, 1)
    END IF
  END IF



  ! ----------------------------------------------------------------------------
  ! --- 3. INTEGRATE: the ACCELERATION regions: areaA ---
  IF ( converged_Pre) THEN
    areaA = 0.0_C_DOUBLE
    count_Acc_regions = 0
  ELSE
    zeroL = zeroR  ! The last region's right-side zero is next region's left-side zero
    CALL integrateAccelerationRegions(m, leftOfMax, zeroL,  tmax,                               & ! INPUTS
                                      areaA, zeroR, count_Acc_regions, converged_Accelerating)    ! OUTPUTS
  
    IF (Cverbose) THEN
      IF ( .NOT.(converged_Accelerating) ) THEN
        CALL DBLEPR(" Accelerating dod not converged by t =", -1, zeroR, 1)
      END IF
      
      CALL DBLEPR("    Aacc area:", -1, areaA, 1)
      CALL DBLEPR("      between:", -1, zeroL, 1)
      CALL DBLEPR("          and:", -1, zeroR, 1)
      CALL INTPR( "      up to m:", -1, m,     1)
    END IF  
  END IF
  
  ! --- WIND THINGS UP ---
  count_Integration_Regions = count_Integration_Regions + count_Acc_regions
  areaT = area0 + area1 + areaA

  IF (Cverbose) THEN
    CALL DBLEPR("* Initial area0: ", -1, area0, 1)
    CALL DBLEPR("* Pre-acc area1: ", -1, area1, 1)
    CALL DBLEPR("*     Acc area!: ", -1, areaA, 1)
    CALL DBLEPR("***       TOTAL: ", -1, areaT, 1)
    CALL INTPR( "   over regions: ", -1, count_Integration_Regions, 1)
  END IF

  ! We have the value of the integral in the PDF/CDF calculation.
  ! So now work out the actualy PDF/CDF
  
  IF (Cpdf) THEN
    funvalueI = areaT/pi 
  ELSE
    funvalueI = -areaT/pi + 0.5E0_C_DOUBLE
  END IF  
  IF (Cverbose) CALL DBLEPR("***    Fun. value:", -1, funvalueI, 1)
  

  CONTAINS
  

    SUBROUTINE advanceM(i, m, mmax, mOld, leftOfMax, flip)
      ! Determine the next value of m, for solving the zeros of the integrand
      
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
    
      IMPLICIT NONE
      
      INTEGER(C_INT), INTENT(IN)    :: mmax, i    ! Maximum value m can take, and current index
      INTEGER(C_INT), INTENT(INOUT) :: m          ! M index (used for calculation and C-binding)
      INTEGER(C_INT), INTENT(OUT)   :: mOld       ! Previous index value
      INTEGER(C_INT), INTENT(INOUT) :: leftOfMax  ! True if on the left side of kmax
      INTEGER(C_INT), INTENT(INOUT) :: flip       ! True if cross from left to right
      
      REAL(KIND=C_DOUBLE)           :: current_y, current_mu, current_phi
    
    
      ! Grab the relevant scalar values for this iteration:
      current_y    = Cy(i)    ! Access y value for index i
      current_mu   = Cmu(i)   ! Access mu value for index i
      current_phi  = Cphi(i)  ! Access phi value for index i
    
      mOld = m
      flip = 0
      
      IF (current_y .GE. current_mu) THEN
        ! Always heading downwards (away from kmax), so easy
        m = m - 1 
      ELSE
        ! We have a maximum (kmax) to consider
        IF (leftOfMax .EQ. 1) THEN
          IF (m == mmax) THEN 
            ! Move to the other side of the maximum
            leftOfMax = 0
            flip = 1
          ELSE
            ! Continue towards the maximum
            m = m + 1 
            ! mOld is already saved before the IF block
            leftOfMax = 1 ! Still on the left side
          END IF
        ELSE
          ! When on the RIGHT of the maximum, can always just reduce m by one
          m = m - 1 
          leftOfMax = 0
        END IF
      END IF
    
    END SUBROUTINE advanceM
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SUBROUTINE integrateFirstRegion(m, leftOfMax, area0, zeroR)
      ! Integrates the initial region of the Fourier integral

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: area0, zeroR
      INTEGER(C_INT), INTENT(IN)          :: m
      INTEGER(C_INT), INTENT(INOUT)       :: leftOfMax
      REAL(KIND=C_DOUBLE)                 :: t_Start_Point, zeroL, zeroBoundL, zeroBoundR
      REAL(KIND=C_DOUBLE)                 :: pi


      pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
      count_Integration_Regions = 1

      ! Find starting point for the first zero
      IF (leftOfMax .EQ. 1) THEN
        t_Start_Point = pi / current_y  
        zeroBoundL = 0
        zeroBoundR = tmax   ! WAS: t_Start_Point * 2.0_C_DOUBLE
      ELSE
        ! Searching to the right of tmax
        t_Start_Point = tmax + pi / current_y  
        zeroBoundL = tmax
        zeroBoundR = t_Start_Point * 2.0_C_DOUBLE
      END IF
    
      ! Find the zero
      zeroL = 0.0_C_DOUBLE
      CALL findExactZeros(i, mfirst, zeroBoundL, zeroBoundR, t_Start_Point, zeroR, leftOfMax)
      CALL evaluateImk(i, zeroR, TMP)

      ! Find the area
      CALL GaussQuadrature(i, zeroL, zeroR, area0)
      
      IF (Cverbose) THEN
        CALL DBLEPR("  *** INITIAL area:", -1, area0, 1 )
        CALL DBLEPR("         between t:", -1, 0.0_C_DOUBLE, 1 )
        CALL DBLEPR("             and t:", -1, zeroR, 1 )
        CALL INTPR( "   using (right) m:", -1, m, 1)
      END IF

      END SUBROUTINE integrateFirstRegion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE integratePreAccRegions(m, mfirst, leftOfMax, zeroL,  tmax,  tStartAcc,                  & ! INPUTS
                                     area1, zeroR, count_PreAcc_regions, converged_Pre)                  ! OUTPUTS
      ! Integration of the region *after* the initial, but *before* acceleration is invoked.
      ! Potentially, everything converges in this step (without ever needing acceleration).
      
      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: area1, zeroR
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax, tStartAcc
      INTEGER(C_INT)                      :: mfirst
      INTEGER(C_INT), INTENT(INOUT)       :: m, leftOfMax
      INTEGER(C_INT), INTENT(OUT)         :: count_PreAcc_regions
      LOGICAL(C_BOOL), INTENT(OUT)        :: converged_Pre

      INTEGER(C_INT)                      :: mOld
      REAL(KIND=C_DOUBLE)                 :: zeroL, zeroBoundL, zeroBoundR
      REAL(KIND=C_DOUBLE)                 :: area1Old, tolerance, sumAOld
      LOGICAL(C_BOOL)                     :: stop_PreAccelerate


      converged_Pre = .FALSE.       ! .TRUE. if it seems the integration has converged
      tolerance = 1.0E-12_C_DOUBLE  ! tolerance for concluding the area is not vhanging
      sumAOld = 0.0_C_DOUBLE        ! Previous estimation of the area
      count_PreAcc_regions = 0      ! Count how many pre-acc regions are evaluated

      IF (mfirst .EQ. -1) THEN
        area1 = 0.0_C_DOUBLE
        mOld = m
        zeroR = zeroL 
    
        IF (Cverbose) CALL DBLEPR("  - No PRE-ACC area for y:", -1, current_y, 1 )
      ELSE
        area1 = 0.0_C_DOUBLE
        area1Old = 10.0_C_DOUBLE
        mOld = m
    
        CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)

        stop_PreAccelerate = .FALSE.
        DO WHILE ( .NOT.(stop_PreAccelerate ) )
          count_PreAcc_regions = count_PreAcc_regions + 1
    
          IF (leftOfMax .EQ. 1) THEN
            zeroBoundL = zeroR
            zeroBoundR = tmax
            IF (flip_To_Other_Side .EQ. 1) THEN 
              leftOfMax = 0
              flip_To_Other_Side = 0
            END IF
          ELSE
            zeroBoundL = tmax
            zeroBoundR = zeroBoundL * 20.0_C_DOUBLE
          END IF
    
          zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
          zeroL = zeroR
    
          ! ONLY NEEDED FOR SMALL P???
          CALL improveKZeroBounds(i, m, leftOfMax, zeroStartPoint, zeroBoundL, zeroBoundR)
          CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zero, leftOfMax)

          zeroR = zero
    
          IF (count_Integration_Regions .GT. 1) THEN
            area1Old = area1
            sumAOld = sumA
          END IF

          CALL GaussQuadrature(i, zeroL, zeroR, sumA)
          area1 = area1 + sumA
          
          count_Integration_Regions = count_Integration_Regions + 1
    
          ! Stop condition for pre-acceleration.
          ! Ensure that we have passed the peak of Im k(t), so that acceleration can be used
          IF ( (count_PreAcc_regions .GE. 2) .AND. (zeroL .GT. tStartAcc) ) THEN
            stop_PreAccelerate = .TRUE.
          END IF
         
          ! Check for convergence
          IF ( (DABS(sumA - sumAOld) .LE. tolerance)  ) THEN
            converged_Pre = .TRUE.
            stop_PreAccelerate = .TRUE.  
          END IF

          mOld = m
          CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)
        END DO 
    
        IF (Cverbose) THEN
          CALL DBLEPR("  *** PRE-ACC area:", -1, sumA, 1 )
          CALL DBLEPR("         between t:", -1, zeroL, 1 )
          CALL DBLEPR("             and t:", -1, zeroR, 1 )
          CALL INTPR( "           using m:", -1, m, 1 )
        END IF
    
      END IF
      
    END SUBROUTINE integratePreAccRegions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE integrateAccelerationRegions(m, leftOfMax, zeroL,  tmax, & ! INPUTS
                                            West, zeroR, its_Acceleration, converged_Accelerating)
      ! Integration of the region *after* the initial and pre-acceleration regions, the tail area.
      ! Potentially, everything converges in pre-acceleration (without needing acceleration).

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroR
      REAL(KIND=C_DOUBLE)                 :: tmax, zeroL
      INTEGER(C_INT), INTENT(INOUT)       :: m, leftOfMax
      INTEGER(C_INT), INTENT(OUT)         :: its_Acceleration

      INTEGER(C_INT)                      :: mOld
      REAL(KIND=C_DOUBLE)                 :: zeroBoundL, zeroBoundR
      REAL(KIND=C_DOUBLE)                 :: West, Wold, Wold2
      LOGICAL(C_BOOL)                     :: keep_Accelerating, converged_Accelerating


      ! Initialisation of acceleration 
      West = 3.0_C_DOUBLE                 ! The current estimate of the tail area
      Wold = 2.0_C_DOUBLE                 ! Wold and Wold2 are the previous two estimates of the tail area,
      Wold2 = 1.0_C_DOUBLE                !   to allow testing for convergence
      its_Acceleration = 0                ! Number of integration regions used
      keep_Accelerating = .TRUE.          ! .TRUE. if we need more integration regions
      converged_Accelerating = .FALSE.    ! .TRUE. if it seems the integration has converged
      
      accMax = 100             ! Maximum number of regions in acceleration; arbitrary
      min_Acc_Regions = 3     ! Minimum number of acceleration regions to use; need at leats three

      ! This will be the very first, left-most value of t used in acceleration
      xvec(1) = zeroL
      
      DO WHILE (keep_Accelerating )
        its_Acceleration = its_Acceleration + 1
    
        mOld = m
        CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)

        zeroBoundL = zeroL
    
        IF (leftOfMax .EQ. 1) THEN
          zeroBoundR     = tmax
          zeroStartPoint = (zeroBoundL + zeroBoundR) / 2.0_C_DOUBLE
        ELSE
          ! Searching to the right of the maximum of Im k(t)
          zeroStartPoint = zeroL + (pi / current_y)
          zeroBoundR = zeroStartPoint * 5.0_C_DOUBLE
          IF (flip_To_Other_Side .EQ. 1) THEN 
            ! Then this is the first time on the right side of the maximum
            leftOfMax = 0
            flip_To_Other_Side = 0
          END IF
        END IF

        ! Find the exact zero
        CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zeroR, leftOfMax)
        IF (its_Acceleration .GE. accMax) THEN
          IF (Cverbose) CALL INTPR("Max acceleration regions reached. Stopping acceleration at m:", -1, m, 1)
          convergence_Acc = 1
          EXIT
        END IF
    
        xvec(its_Acceleration + 1) = zeroR
        
        CALL GaussQuadrature(i, zeroL, zeroR, psi)           ! psi: area of the latest region
        count_Integration_Regions = count_Integration_Regions + 1

        wvec(its_Acceleration) = psi 
          ! wvec contains the sequence of integration areas, starting with the first and up to the limit
    
        ! Update past estimates
        Wold2 = Wold
        Wold  = West    
        
        ! Estimate the tail area using acceleration
        CALL accelerate(xvec, wvec, its_Acceleration, Mmatrix, Nmatrix, West)

        ! Check for convergence
        relerr = ( DABS(West - Wold) + DABS(West - Wold2)) / (DABS(West) + epsilon)
        IF (Cverbose) THEN
          CALL DBLEPR("  *** This acceleration area:", -1, psi, 1)
          CALL DBLEPR("                   between t:", -1, zeroL, 1)
          CALL DBLEPR("                       and t:", -1, zeroR, 1)
          CALL INTPR( "                     using m:", -1, m, 1)
          CALL DBLEPR("                with rel err:", -1, relerr, 1)
        END IF
    
        ! Declare convergence_Acc of we have sufficient regions, and relerr estimate is small
        IF ( (its_Acceleration .GE. min_Acc_Regions) .AND. &
             (relerr .LT. aimrerr) ) THEN
          keep_Accelerating = .FALSE.
          converged_Accelerating = .TRUE.
          convergence_Acc = 1
        END IF
        
        IF (its_Acceleration .GT. accMax) THEN
             keep_Accelerating = .FALSE.
             converged_Accelerating = .FALSE.
        END IF
!IF (m .LT. -13) STOP

        zeroL = zeroR ! Next left zero, is previous region's right zero

      END DO  
    
    END SUBROUTINE integrateAccelerationRegions

END SUBROUTINE ComputeTwIntegral

