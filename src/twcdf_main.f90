SUBROUTINE twcdf_main(N, p, phi, y, mu, verbose, funvalue, exitstatus, relerr, Int_Regions)
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)  :: N, verbose
  REAL(C_DOUBLE), INTENT(IN)  :: p
  REAL(C_DOUBLE), INTENT(IN)  :: phi(N), y(N), mu(N)
  REAL(C_DOUBLE), INTENT(OUT) :: funvalue(N)
  INTEGER(C_INT), INTENT(OUT) :: exitstatus
  REAL(C_DOUBLE), INTENT(OUT) :: relerr
  INTEGER(C_INT), INTENT(OUT) :: Int_Regions
  ! --- EXPLICIT INTERFACES FOR INTERNAL CALLS ---
  INTERFACE
    SUBROUTINE DFcompute(i, funvalueI, exitstatus, relerr, Int_Regions)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: i
      REAL(KIND=C_DOUBLE), INTENT(INOUT)        :: funvalueI
      INTEGER, INTENT(OUT)                      :: exitstatus
      REAL(KIND=C_DOUBLE), INTENT(OUT)          :: relerr
      INTEGER, INTENT(OUT)                      :: Int_Regions
    END SUBROUTINE DFcompute
  END INTERFACE
    ! -----------------------------------------------

  ! --- Local variables ---
  INTEGER               :: i
  REAL(KIND=C_DOUBLE)   :: funvalueTMP
  

  ! --- Initialization ---
  Cp = p
  Cy = y
  Cmu = mu
  Cphi = phi
  CN = N
  IF (verbose .EQ. 1) THEN
    Cverbose = .TRUE.
  ELSE
    Cverbose = .FALSE.
  END IF

  exitstatus = 1
  relerr = 0.0_C_DOUBLE


  ! --- Determine case: psmall = TRUE means 1 < p < 2 ---
  CpSmall = .FALSE.
  IF ( (p > 1.0_C_DOUBLE) .AND. (p < 2.0_C_DOUBLE) ) CpSmall = .TRUE.

  ! --- Loop over N values ---
  DO i = 1, N
    CALL DFcompute(i, funvalueTMP, exitstatus, relerr, Int_Regions)
    funvalue(i) = funvalueTMP
  END DO

END SUBROUTINE twcdf_main
