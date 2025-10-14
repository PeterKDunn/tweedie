
REAL(KIND=C_DOUBLE) FUNCTION findImkM(i, t, m) RESULT(ImkM_out) BIND(C, NAME='findImkM')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  
  ! --- Interface Declaration ---
  ! We declare the interface for the subroutine we must call (findImk)
  INTERFACE
    SUBROUTINE findImk(i, t_in, ImkM_out)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t_in
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: ImkM_out
    END SUBROUTINE findImk
  END INTERFACE
  
  REAL(KIND=C_DOUBLE), INTENT(IN) :: t
  REAL(KIND=C_DOUBLE)             :: pi, Imk_val
  INTEGER(C_INT), INTENT(IN)      :: m, i

  pi = 4.0D0 * DATAN(1.0D0)

  CALL findImk(i, t, Imk_val)
  
  ImkM_out = Imk_val - REAL(m, KIND=C_DOUBLE) * pi

END FUNCTION findImkM
