
SUBROUTINE findImkdZero(i, t, f, df) 
  USE tweedie_params_mod, ONLY: Cphi, Cmu, Cy
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
  REAL(KIND=C_DOUBLE)  :: current_y, current_mu, current_phi
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  CALL findImkd(i, t, Imkd)
  CALL findImkdd(i, t, Imkdd)

  f  = Imkd
  df = Imkdd

END SUBROUTINE findImkdZero
