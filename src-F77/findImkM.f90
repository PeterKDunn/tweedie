
SUBROUTINE findImkM(t, f, df)
  
  USE tweedie_params_mod
  
  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN) :: t
  DOUBLE PRECISION, INTENT(OUT) :: f, df
  DOUBLE PRECISION :: pi
  DOUBLE PRECISION :: Imk, Imdk
  
  ! --- Executable Statements ---

  CALL findImk(t, Imk)
  CALL findImkd(t, Imdk)
  
  f  = Imk - DBLE(m) * PI
  df = Imdk
  
END SUBROUTINE findImkM
