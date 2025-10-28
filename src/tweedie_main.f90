SUBROUTINE twcomputation_main(N, p, phi, y, mu, verbose, pdf, funvalue, exitstatus, relerr, Int_Regions)
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)  :: N, verbose, pdf
  REAL(C_DOUBLE), INTENT(IN)  :: p
  REAL(C_DOUBLE), INTENT(IN)  :: phi(N), y(N), mu(N)
  REAL(C_DOUBLE), INTENT(OUT) :: funvalue(N)
  INTEGER(C_INT), INTENT(OUT) :: exitstatus
  REAL(C_DOUBLE), INTENT(OUT) :: relerr
  INTEGER(C_INT), INTENT(OUT) :: Int_Regions(N)
  ! --- EXPLICIT INTERFACES FOR INTERNAL CALLS ---
  INTERFACE
    SUBROUTINE IntegralCompute(i, funvalueI, exitstatus, relerr, Int_Regions)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      USE tweedie_params_mod
      
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: i
      REAL(KIND=C_DOUBLE), INTENT(INOUT)        :: funvalueI
      INTEGER, INTENT(OUT)                      :: exitstatus
      REAL(KIND=C_DOUBLE), INTENT(OUT)          :: relerr
      INTEGER, INTENT(OUT)                      :: Int_Regions
    END SUBROUTINE IntegralCompute
  END INTERFACE
    ! -----------------------------------------------

  ! --- Local variables ---
  INTEGER               :: i, Int_RegionsTMP
  REAL(KIND=C_DOUBLE)   :: funvalueTMP
  

  ! --- Initialization ---
  Cp = p
  Cy = y
  Cmu = mu
  Cphi = phi
  CN = N
  IF (pdf .EQ. 0) THEN
    ! Coputing the DF
    Cpdf = .FALSE.
  ELSE
    ! Computing the PDF
    Cpdf = .TRUE.
  END IF
  IF (verbose .EQ. 1) THEN
    ! Verbose feedback
    Cverbose = .TRUE.
  ELSE
    ! Minimal feedback
    Cverbose = .FALSE.
  END IF

  exitstatus = 1
  relerr = 0.0_C_DOUBLE


  ! --- Determine case: psmall = TRUE means 1 < p < 2 ---
  CpSmall = .FALSE.
  IF ( (p > 1.0_C_DOUBLE) .AND. (p < 2.0_C_DOUBLE) ) CpSmall = .TRUE.

  ! --- Loop over N values ---
  DO i = 1, N
    CALL IntegralCompute(i, funvalueTMP, exitstatus, relerr, Int_RegionsTMP)
    funvalue(i) = funvalueTMP
    Int_Regions(i) = Int_RegionsTMP
  END DO

END SUBROUTINE twcomputation_main
