SUBROUTINE twpdf_main(N, p, phi, y, mu, exact, funvalue, exitstatus, relerr, its)
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER(C_INT), INTENT(IN) :: N
  REAL(C_DOUBLE), INTENT(IN) :: p
  INTEGER(C_INT), INTENT(IN) :: exact
  REAL(C_DOUBLE), INTENT(IN) :: phi(N), y(N), mu(N)
  REAL(C_DOUBLE), INTENT(OUT):: funvalue(N)
  INTEGER(C_INT), INTENT(OUT):: exitstatus, its
  REAL(C_DOUBLE), INTENT(OUT):: relerr
  

  ! --- EXPLICIT INTERFACES FOR INTERNAL CALLS ---
  INTERFACE
    SUBROUTINE PDFbigp(i, exact, funvalueI, exitstatus, relerr, verbose)
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: i, exact
      REAL(KIND=8), INTENT(INOUT)               :: funvalueI
      INTEGER, INTENT(OUT)                      :: exitstatus
      REAL(KIND=8), INTENT(OUT)                 :: relerr
      INTEGER, INTENT(IN)                       :: verbose
    END SUBROUTINE PDFbigp

    SUBROUTINE PDFsmallp(i, exact, funvalueI, exitstatus, relerr, verbose)
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: i, exact
      REAL(KIND=8), INTENT(INOUT)               :: funvalueI
      INTEGER, INTENT(OUT)                      :: exitstatus
      REAL(KIND=8), INTENT(OUT)                 :: relerr
      INTEGER, INTENT(IN)                       :: verbose
    END SUBROUTINE PDFsmallp
  END INTERFACE


  INTEGER         :: i
  REAL(C_DOUBLE)  :: aimrerr, funvalueTMP
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
  aimrerr = 1.0E-10_C_DOUBLE

  ! --- Determine case: psmall = TRUE means 1 < p < 2 ---
  CpSmall = .FALSE.
  IF ( (p > 1.0_C_DOUBLE) .AND. (p < 2.0_C_DOUBLE) ) CpSmall = .TRUE.

  ! --- Loop over N ---
  DO i = 1, N
    IF (CpSmall) THEN
      CALL PDFsmallp(i, exact, funvalueTMP, exitstatus, relerr, verbose)
    ELSE
      CALL PDFbigp(i, exact, funvalueTMP, exitstatus, relerr, verbose)
    END IF
    funvalue(i) = funvalueTMP
  END DO

END SUBROUTINE twpdf_main

