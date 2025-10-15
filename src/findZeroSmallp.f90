
SUBROUTINE findZeroSmallp(i, t, f, df)
  USE tweedie_params_mod

  IMPLICIT NONE
  
  ! Variables f and df are the function value and derivative output
  REAL(KIND=8), INTENT(IN)    :: t
  INTEGER, INTENT(IN)         :: i
  REAL(KIND=8), INTENT(OUT)   :: f, df
  
  ! --- Local Variables ---
  REAL(KIND=8) :: Imk, Imdk
  REAL(KIND=8) :: current_y, current_mu, current_phi
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  CALL findImk(i, t, Imk)
  CALL findImkd(i, t, Imdk)
  
  ! --- Core Calculation ---
  f  = DSIN(Imk + (t * current_y) ) - DSIN( Imk)
  df = (Imdk + current_y) * DCOS(Imk + (t * current_y) ) - &
       Imdk * DCOS(Imk)
  
END SUBROUTINE findZeroSmallp
