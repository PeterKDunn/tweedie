
SUBROUTINE findLambda(i, lambda)
  USE tweedie_params_mod

  IMPLICIT NONE
  
  INTEGER, INTENT(IN)        :: i
  REAL(KIND=8), INTENT(OUT)  :: lambda ! The output value
  
  REAL(KIND=8)                :: current_mu, current_phi

  ! Grab the relevant scalar values for this iteration:
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  

  lambda = 0.0d0
  IF (CpSmall) THEN
    ! The calculation for lambda (used in P(Y=0) = exp(-lambda))
    lambda = (current_mu ** (2.0d00 - Cp) ) / &
             (current_phi * (2.0d00 - Cp) )
    ! NOTE: No negative sign in front
  END IF
  
END SUBROUTINE findLambda
