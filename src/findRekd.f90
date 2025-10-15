
SUBROUTINE findRekd(i, t, Redk)
  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu

  IMPLICIT NONE
  
  ! Arguments (Inputs/Outputs)
  REAL(KIND=8), INTENT(IN)    :: t              ! Input parameter t
  INTEGER, INTENT(IN)         :: i
  REAL(KIND=8), INTENT(OUT)   :: Redk           ! Output result (Derivative of Real part of k(t))
  
  REAL(KIND=8)            :: omega, pindex
  REAL(KIND=8)            :: current_mu, current_phi

  ! Grab the relevant scalar values for this iteration:
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  pindex = 1.0d00 / (1.0d00 - Cp)
  
  omega = ATAN( ( (1.0d00 - Cp) * t * current_phi) / &
                (current_mu ** (1.0d00 - Cp) ) )
  
  Redk = current_mu * &
         SIN( omega * pindex ) / &
         (COS(omega)**pindex)
  
END SUBROUTINE findRekd
