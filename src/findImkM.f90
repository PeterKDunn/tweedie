
REAL(KIND=8) FUNCTION findImkM(i, t, m) RESULT(ImkM_out) 
  USE tweedie_params_mod

  IMPLICIT NONE

  
  ! --- Interface Declaration ---
  ! We declare the interface for the subroutine we must call (findImk)
  INTERFACE
    SUBROUTINE findImk(i, t_in, ImkM_out)

      IMPLICIT NONE
      INTEGER, INTENT(IN)        :: i
      REAL(KIND=8), INTENT(IN)   :: t_in
      REAL(KIND=8), INTENT(OUT)  :: ImkM_out
    END SUBROUTINE findImk
  END INTERFACE
  
  REAL(KIND=8), INTENT(IN) :: t
  REAL(KIND=8)             :: pi, Imk_val
  INTEGER, INTENT(IN)      :: m, i

  pi = 4.0D0 * DATAN(1.0D0)

  CALL findImk(i, t, Imk_val)
  
  ImkM_out = Imk_val - REAL(m, KIND=8) * pi

END FUNCTION findImkM
