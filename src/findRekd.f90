SUBROUTINE findRekd(i, t, Redk)
  ! Find the value of Re k'(t)
  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN)    :: t
  INTEGER(C_INT), INTENT(IN)         :: i
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: Redk
  
  REAL(KIND=C_DOUBLE)       :: omega, pindex
  REAL(KIND=C_DOUBLE)       :: current_mu, current_phi


  ! Grab the relevant scalar values for this iteration:
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  
  pindex = 1.0E0_C_DOUBLE / (1.0E0_C_DOUBLE - Cp)
  omega = DATAN( ( (1.0E0_C_DOUBLE - Cp) * t * current_phi) / &
                 (current_mu ** (1.0E0_C_DOUBLE - Cp) ) )
  
  Redk = current_mu * &
         DSIN( omega * pindex ) / &
         (DCOS(omega)**pindex)
  
END SUBROUTINE findRekd
