SUBROUTINE twpdf(N, p, phi, y, mu, exact, funvalue, exitstatus, relerr, its) BIND(C, name="twpdf")
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE tweedie_params_mod

  IMPLICIT NONE
  INTEGER(C_INT), INTENT(IN)  :: N
  REAL(C_DOUBLE), INTENT(IN)  :: p, phi(N), y(N), mu(N)
  REAL(C_DOUBLE), INTENT(OUT) :: funvalue(N), relerr
  INTEGER(C_INT), INTENT(IN)  :: exact
  INTEGER(C_INT), INTENT(OUT) :: exitstatus, its

  ! Call internal Fortran routine
  CALL twpdf_main(N, p, phi, y, mu, exact, funvalue, exitstatus, relerr, its)
  
END SUBROUTINE twpdf
