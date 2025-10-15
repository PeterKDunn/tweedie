SUBROUTINE findRekd(i, t, Redk)
  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Arguments (Inputs/Outputs)
  REAL(KIND=C_DOUBLE), INTENT(IN)    :: t              ! Input parameter t
  INTEGER(C_INT), INTENT(IN)         :: i
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: Redk           ! Output result (Derivative of Real part of k(t))
  
  REAL(KIND=C_DOUBLE)            :: omega, pindex
  REAL(KIND=C_DOUBLE)            :: current_mu, current_phi

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
