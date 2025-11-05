SUBROUTINE evaluateLambda(i, lambda)
  ! Find lambda, such that P(Y = 0) = exp( -lambda ) when 1 < p < 2 
  
  USE tweedie_params_mod, ONLY: Cmu, Cphi, Cp, CpSmall
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: lambda 
  
  REAL(KIND=C_DOUBLE)               :: current_mu, current_phi


  ! Grab the relevant scalar values for this iteration:
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  

  lambda = 0.0E0_C_DOUBLE
  IF (CpSmall) THEN
    ! The calculation for lambda (used in P(Y=0) = exp(-lambda))
    lambda = (current_mu ** (2.0E0_C_DOUBLE - Cp) ) / &
             (current_phi * (2.0E0_C_DOUBLE - Cp) )
    ! NOTE: No negative sign in front
  END IF
  
END SUBROUTINE evaluateLambda
