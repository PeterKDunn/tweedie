SUBROUTINE twcomputation(N, p, phi, y, mu, verbose, pdf, tail, funvalue, exitstatus, relerr, its) BIND(C, name="twcomputation")
  ! Called by R; calls FORTRAN
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  INTEGER(C_INT), INTENT(IN)        :: N, verbose, pdf, tail
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: p, phi(N), y(N), mu(N)
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalue(N), relerr(n)
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus(N), its(N)


  ! Call internal Fortran routine
  CALL twcomputation_loop(N, p, phi, y, mu, verbose, pdf, tail,   & ! INPUTS
                          funvalue, exitstatus, relerr, its)        ! OUTPUTS

END SUBROUTINE twcomputation
