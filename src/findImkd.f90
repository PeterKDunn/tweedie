
SUBROUTINE findImkd(i, t, Imdk) 

  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu, Cy

  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN)        :: i
  REAL(KIND=8), INTENT(IN)   :: t
  REAL(KIND=8), INTENT(OUT)  :: Imdk  ! The result of the calculation
  
  ! Local Variables
  REAL(KIND=8) :: omega, pindex
  REAL(KIND=8) :: current_y, current_mu, current_phi

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  ! --- Calculation (Using F90 constants) ---
  
  pindex = 1.0_8 / (1.0_8 - Cp)
  
  ! Use DATAN for double precision arc tangent
  omega = DATAN( ( (1.0_8 - Cp) * t * current_phi) / (current_mu ** (1.0_8 - Cp) ) )
  
  ! Final calculation
  Imdk = current_mu * (DCOS(omega * pindex) / (DCOS(omega) ** pindex)) - current_y

END SUBROUTINE findImkd
