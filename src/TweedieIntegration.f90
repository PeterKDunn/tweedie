
SUBROUTINE TweedieIntegration(i, funvalueI, exitstatus, relerr, count_Integration_Regions) 
  ! Compute the value of the integrals in the Fourier-inversion expressions for the PDF and CDF

  USE Integrands_MOD, ONLY: Integrands
  USE tweedie_params_mod
  USE TweedieIntZones
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE rprintf_mod
  USE Calcs_Imag
  USE Calcs_Real
  USE Calcs_K
  
  IMPLICIT NONE
  
  INTEGER(C_INT), INTENT(IN)        :: i              ! Observation index
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalueI      ! Computed result
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr         ! Estimate of relative error
  INTEGER(C_INT), INTENT(OUT)       :: count_Integration_Regions ! Num int regions


  ! Local Variables: All local variables defined here
  INTEGER(C_INT)        :: mmax, mfirst, mOld, accMax
  INTEGER(C_INT)        :: count_PreAcc_Regions, count_Acc_Regions
  INTEGER(C_INT)        :: m, min_Acc_Regions
  LOGICAL(C_BOOL)       :: converged_Accelerating, converged_Pre, convergence_Acc
  REAL(KIND=C_DOUBLE)   :: kmax, tmax, aimrerr
  REAL(KIND=C_DOUBLE)   :: epsilon, areaT, pi, zero
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR, area0, area1, areaA
  REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE), ALLOCATABLE   :: Mmatrix(:, :), Nmatrix(:, :), xvec(:), wvec(:)
  REAL(KIND=C_DOUBLE)   :: leftPreAccZero, leftAccZero
  REAL(KIND=C_DOUBLE)   :: zeroStartPoint
  LOGICAL(C_BOOL)       :: left_Of_Max, flip_To_Other_Side
  

  INTERFACE

    SUBROUTINE GaussQuadrature(i, a, b, area)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(KIND=C_DOUBLE), INTENT(OUT)      :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)       :: a, b
      INTEGER(C_INT), INTENT(IN)            :: i
    END SUBROUTINE GaussQuadrature


    SUBROUTINE accelerate(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)                 :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: xvec(:), wvec(:), Mmatrix(:, :), Nmatrix(:, :)
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: West
    END SUBROUTINE accelerate

  END INTERFACE


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! ALLOCATE these arrays onto the HEAP
  ALLOCATE(Mmatrix(2, 502))
  ALLOCATE(Nmatrix(2, 502))
  ALLOCATE(xvec(502))
  ALLOCATE(wvec(502))
  
  IF (Cverbose) THEN
    CALL DBLEPR("*** Computing for p =", -1, Cp, 1)
    CALL DBLEPR("*** Computing for y =", -1, current_y, 1)
    CALL DBLEPR("*** Computing for mu =", -1, current_mu, 1)
    CALL DBLEPR("*** Computing for phi =", -1, current_phi, 1)
  END IF


  ! --- Initialization ---
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  aimrerr = 1.0E-12_C_DOUBLE ! 1.0E-12_C_DOUBLE
  mOld = 0
  m = 0
  exitstatus = 0
  relerr = 1.0_C_DOUBLE
  convergence_Acc = .FALSE.
  epsilon = 1.0E-12_C_DOUBLE
  Mmatrix = 0.0_C_DOUBLE
  Nmatrix = 0.0_C_DOUBLE
  xvec = 0.0_C_DOUBLE
  wvec = 0.0_C_DOUBLE
  mmax = 0_C_INT
  accMax = 500_C_INT                    ! Max acceleration regions
  min_Acc_Regions = 3_C_INT             ! Min preacceleration regions
  count_Integration_Regions = 0_C_INT   ! Counter for number of integration regions
  zeroStartPoint = 0.0_C_DOUBLE
  flip_To_Other_Side = .FALSE.
  zeroR = 0.0_C_DOUBLE
  zeroL = 0.0_C_DOUBLE

  ! --- Integration initialization ---
  area0 = 0.0_C_DOUBLE
  area1 = 0.0_C_DOUBLE
  zero  = 0.0_C_DOUBLE

  ! --- Find kmax, tmax, mmax ---
  CALL findKmax(i, kmax, tmax, mmax, mfirst, left_Of_Max)
  IF (Cverbose) THEN
    CALL DBLEPR("  -            kmax:", -1, kmax, 1 )
    CALL DBLEPR("  -            tmax:", -1, tmax, 1 )
    CALL INTPR( "  -            mmax:", -1, mmax, 1 )
    CALL INTPR( "  - first zero at m:", -1, mfirst, 1 )
  END IF
  
  m = mfirst
  
  ! INTEGRATION
  ! Three integration zones:
  !   1. The *initial* area, which is not between zeros of Im{k(t)}: area0
  !   2. The initial area *before* Sidi acceleration is invoked: area1
  !      (For instance, wait until Im{k(t)} is on the downturn.)
  !   3. The area thereafter, upon which Sidi acceleration is
  !      applied; the area returned by acceleration is areaA


  ! ----------------------------------------------------------------------------
  ! --- 1. INTEGRATE FIRST (sometimes non-standard) REGION: area0 ---
!WRITE(*,*) "Initial area"
  CALL integrateFirstRegion(i, mfirst, left_Of_Max, tmax, &  ! INPUTS
                            area0, zeroR)                             ! OUTPUTS

  IF (Cverbose) THEN
    CALL DBLEPR("Initial region area:", -1, area0, 1)
    CALL DBLEPR("      between 0 and:", -1, zeroR, 1)
    CALL INTPR( "      using right m:", -1, m, 1)
  END IF




  ! ----------------------------------------------------------------------------
  ! --- 2. INTEGRATE: the PRE-ACCELERATION regions: area1 ---
!WRITE(*,*) "Starting pre-acceleration"
  zeroL = zeroR  ! The last region's right-side zero is next region's left-side zero
  leftPreAccZero = zeroL
  CALL integratePreAccRegions(i, m, mfirst, left_Of_Max, zeroL,  tmax, mmax, flip_To_Other_Side,    & ! INPUTS
                              area1, zeroR, count_PreAcc_regions,  converged_Pre)     ! OUTPUTS
  count_Integration_Regions = count_Integration_Regions + count_PreAcc_regions

  IF (Cverbose) THEN
    CALL DBLEPR(" Pre-acc area:", -1, area1, 1)
    CALL DBLEPR("      between:", -1, leftPreAccZero, 1)
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
!WRITE(*,*) "Starting acceleration"
  IF ( converged_Pre) THEN
    areaA = 0.0_C_DOUBLE
    count_Acc_Regions = 0_C_INT
  ELSE
    leftAccZero = zeroL
    zeroL = zeroR  ! The last region's right-side zero is next region's left-side zero

    CALL integrateAccelerationRegions(i, m, left_Of_Max, zeroL,  tmax, accMax,                  & ! INPUTS
                                          aimrerr, epsilon, xvec, wvec, Mmatrix, Nmatrix, & 
                                          mmax, min_Acc_Regions, exitstatus, flip_To_Other_Side, &
                                          convergence_Acc, &
                                          zeroR, converged_Accelerating, &
                                          relerr)
    IF (Cverbose) THEN
      IF ( .NOT.(converged_Accelerating) ) THEN
        CALL DBLEPR(" Accelerating dod not converged by t =", -1, zeroR, 1)
      END IF
      
      CALL DBLEPR("     Acc area:", -1, areaA, 1)
      CALL DBLEPR("      between:", -1, leftAccZero, 1)
      CALL DBLEPR("          and:", -1, zeroR, 1)
      CALL INTPR( "      up to m:", -1, m,     1)
    END IF  
  END IF
  

  ! --- WIND THINGS UP ---
  count_Integration_Regions = 1_C_INT  +              &   ! Initial zone has one integration region
                              count_PreAcc_Regions +  &   ! Pre-acc regions
                              count_Acc_Regions           ! Acc regions
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
    funvalueI =  0.5_C_DOUBLE - areaT/pi
  END IF  
  IF (Cverbose) CALL DBLEPR("***    Fun. value:", -1, funvalueI, 1)
  
  ! Deallocate the arrays
  DEALLOCATE(Mmatrix, Nmatrix, xvec, wvec)
  
END SUBROUTINE TweedieIntegration

