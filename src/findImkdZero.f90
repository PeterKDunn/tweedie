SUBROUTINE evaluateImkdZero(i, t, f, df) 
  ! Evaluate Im k'(t)  and  Im k''(t)  for solving for Kmax (i.e., Im k'(t) = 0)
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  INTERFACE
    SUBROUTINE evaluateImkd(i, t, Imkd)
      ! Evaluate Im k'(t)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkd
    END SUBROUTINE evaluateImkd

    
    SUBROUTINE evaluateImkdd(i, t, Imddk)
      ! Evaluate Im k''(t)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imddk
    END SUBROUTINE evaluateImkdd
  END INTERFACE
  

  INTEGER(C_INT), INTENT(IN)          :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: f, df

  REAL(KIND=C_DOUBLE)  :: Imkd, Imkdd


  CALL evaluateImkd( i, t, Imkd)
  CALL evaluateImkdd(i, t, Imkdd)

  f  = Imkd
  df = Imkdd

END SUBROUTINE evaluateImkdZero
