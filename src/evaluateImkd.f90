SUBROUTINE evaluateImkd(i, t, Imkd) 
  ! Evaluate Im k'(t)
  
  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu, Cy
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkd  ! The result of the calculation
  
  REAL(KIND=C_DOUBLE) :: omega, pindex
  REAL(KIND=C_DOUBLE) :: current_y, current_mu, current_phi


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i


  pindex = 1.0_C_DOUBLE / (1.0_C_DOUBLE - Cp)
  omega = DATAN( ( (1.0_C_DOUBLE - Cp) * t * current_phi) / (current_mu ** (1.0_C_DOUBLE - Cp) ) )
  
  Imkd = current_mu * (DCOS(omega * pindex) / (DCOS(omega) ** pindex)) - current_y

END SUBROUTINE evaluateImkd
