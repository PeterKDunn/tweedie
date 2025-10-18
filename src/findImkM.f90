SUBROUTINE findImkM(i, t, f, df, m)
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE findImk(i, t_in, ImkM_out)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t_in
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: ImkM_out
    END SUBROUTINE findImk
  END INTERFACE
  
  
  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
  INTEGER(C_INT), INTENT(IN)        :: m

  REAL(KIND=8)                      :: pi, Imk_val

  pi = 4.0D0 * DATAN(1.0D0)

  CALL findImk(i, t, Imk_val)
  f = Imk_val - REAL(m, KIND=C_DOUBLE) * pi
  CALL findImkd(i, t, df)
  
END SUBROUTINE findImkM
