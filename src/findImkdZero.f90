
SUBROUTINE findImkdZero(i, t, f, df) 
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  INTERFACE
    SUBROUTINE findImkd(i, t, Imkd)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkd
    END SUBROUTINE findImkd

    
    SUBROUTINE findImkdd(i, t, Imddk)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imddk
    END SUBROUTINE findImkdd
  END INTERFACE
  

  REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: f, df
  INTEGER(C_INT), INTENT(IN)          :: i

  REAL(KIND=C_DOUBLE)  :: Imkd, Imkdd


  CALL findImkd(i, t, Imkd)
  CALL findImkdd(i, t, Imkdd)

  f  = Imkd
  df = Imkdd

END SUBROUTINE findImkdZero
