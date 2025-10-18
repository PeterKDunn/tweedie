SUBROUTINE findZeroSmallp(i, t, f, df)
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Variables f and df are the function value and derivative output
  REAL(KIND=C_DOUBLE), INTENT(IN)    :: t
  INTEGER(C_INT), INTENT(IN)         :: i
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: f, df
  
  ! --- Local Variables ---
  REAL(KIND=C_DOUBLE) :: Imk, Imkd
  REAL(KIND=C_DOUBLE) :: current_y, current_mu, current_phi
  
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  CALL findImk(i, t, Imk)
  CALL findImkd(i, t, Imkd)
  
  f  = DSIN(Imk + (t * current_y) ) - DSIN( Imk)
  df = (Imkd + current_y) * DCOS(Imk + (t * current_y) ) - &
        Imkd * DCOS(Imk)
  
END SUBROUTINE findZeroSmallp
