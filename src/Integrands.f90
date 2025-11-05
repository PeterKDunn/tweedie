MODULE Integrands_MOD
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
  IMPLICIT NONE
  
CONTAINS

  FUNCTION Integrands(i, t) RESULT(integrand_result)  BIND(C, NAME='Integrands')
    ! Function to return the integrand values

    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
    IMPLICIT NONE
    
    INTEGER(C_INT), INTENT(IN)        :: i
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t                ! The internal variable for integration

    REAL(KIND=C_DOUBLE)               :: integrand_result ! The result of the function
    REAL(KIND=C_DOUBLE)               :: current_y, current_mu
    REAL(KIND=C_DOUBLE)               :: Imk, Rek, lambda
    
    
    INTERFACE
    
      SUBROUTINE evaluateImk(i, t, Imk)
        ! Find Im k(t)
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
        IMPLICIT NONE

        INTEGER(C_INT), INTENT(IN)          :: i
        REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
        REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Imk
      END SUBROUTINE evaluateImk

      
      SUBROUTINE evaluateRek(i, t, Rek)
        ! Find Re k(t)
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
        IMPLICIT NONE

        INTEGER(C_INT), INTENT(IN)          :: i
        REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
        REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Rek
      END SUBROUTINE evaluateRek
      
      
      SUBROUTINE evaluateLambda(i, lambda)
       ! Find lambda, such that P(Y = 0) = exp( -lambda ) when 1 < p < 2 
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
        
        IMPLICIT NONE

        INTEGER(C_INT), INTENT(IN)        :: i
        REAL(KIND=C_DOUBLE), INTENT(OUT)  :: lambda
      END SUBROUTINE evaluateLambda
      
    END INTERFACE
      
  
    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)
    current_mu   = Cmu(i)


    ! Check for when t = 0 (t is very close to zero)
    IF (ABS(t) .LT. 1.0E-14_C_DOUBLE) THEN
      ! This should ideally be handled by the integrator (limits), 
      ! but returning the analytic limit is safest.
      integrand_result = current_mu - current_y
  
      RETURN
    ELSE
      CALL evaluateRek(i, t, Rek)
      CALL evaluateImk(i, t, Imk)
      
      IF (Cpdf) THEN
        IF (CpSmall) THEN
          CALL evaluateLambda(i, lambda)
          integrand_result = DEXP( Rek ) * DCOS( Imk ) - DEXP( -lambda ) * DCOS(t * current_y )
        ELSE
          integrand_result = DEXP( Rek ) * DCOS( Imk )
        END IF
      ELSE
        integrand_result = DEXP( Rek ) * DSIN( Imk ) / t
      END IF
    END IF
    
  END FUNCTION Integrands

END MODULE Integrands_MOD
