SUBROUTINE twcdf_main(N, p, phi, y, mu, funvalue, exitstatus, relerr, its)
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)  :: N
  REAL(C_DOUBLE), INTENT(IN)  :: p
  REAL(C_DOUBLE), INTENT(IN)  :: phi(N), y(N), mu(N)
  REAL(C_DOUBLE), INTENT(OUT) :: funvalue(N)
  INTEGER(C_INT), INTENT(OUT) :: exitstatus
  REAL(C_DOUBLE), INTENT(OUT) :: relerr
  INTEGER(C_INT), INTENT(OUT) :: its

  ! --- EXPLICIT INTERFACES FOR INTERNAL CALLS ---
  INTERFACE
    SUBROUTINE DFbigp(i, funvalue, exitstatus, relerr, verbose)
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: i
      REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: funvalue
      INTEGER, INTENT(OUT)                      :: exitstatus
      REAL(KIND=8), INTENT(OUT)                 :: relerr
      INTEGER, INTENT(IN)                       :: verbose
    END SUBROUTINE DFbigp

    SUBROUTINE DFsmallp(i, funvalue, exitstatus, relerr, verbose)
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: i
      REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: funvalue
      INTEGER, INTENT(OUT)                      :: exitstatus
      REAL(KIND=8), INTENT(OUT)                 :: relerr
      INTEGER, INTENT(IN)                       :: verbose
    END SUBROUTINE DFsmallp
  END INTERFACE
    ! -----------------------------------------------

  ! --- Local variables ---
  INTEGER         :: i
  INTEGER         :: verbose

  ! --- Initialization ---
  Cp = p
  Cy = y
  Cmu = mu
  Cphi = phi
  CN = N
  verbose = 1
  exitstatus = 1
  relerr = 0.0_C_DOUBLE
  its = 0

  ! --- Determine case: psmall = TRUE means 1 < p < 2 ---
  CpSmall = .FALSE.
  IF ( (p > 1.0_C_DOUBLE) .AND. (p < 2.0_C_DOUBLE) ) CpSmall = .TRUE.

  ! --- Loop over N values ---
  DO i = 1, N
    IF (CpSmall) THEN
      CALL DFsmallp(i, funvalue, exitstatus, relerr, verbose)
    ELSE
      CALL DFbigp(i, funvalue, exitstatus, relerr, verbose)
    END IF
  END DO

END SUBROUTINE twcdf_main
