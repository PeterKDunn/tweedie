MODULE DFintegrand_MOD
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  IMPLICIT NONE
  
CONTAINS
  FUNCTION DFintegrand(i, t) RESULT(integrand_result)  BIND(C, NAME='dfintegrand')
    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
    IMPLICIT NONE
    
    ! Arguments (Inputs)
    INTEGER(C_INT), INTENT(IN)        :: i
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t                ! The internal variable for integration
    REAL(KIND=C_DOUBLE)               :: integrand_result ! The result of the function
    REAL(KIND=C_DOUBLE)               :: current_y, current_mu
    REAL(KIND=C_DOUBLE)               :: Imk, Rek
    
    
    INTERFACE
      SUBROUTINE findImk(i, t, Imk)
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
        IMPLICIT NONE
        INTEGER(C_INT), INTENT(IN)          :: i
        REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
        REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Imk
      END SUBROUTINE findImk

      
      SUBROUTINE findRek(i, t, Rek)
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
        IMPLICIT NONE
        INTEGER(C_INT), INTENT(IN)          :: i
        REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
        REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Rek
      END SUBROUTINE findRek
    END INTERFACE
      
  
    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)
    current_mu   = Cmu(i)

    ! Check for when t = 0 (t is very close to zero)
    IF (ABS(t) < 1.0E-14_C_DOUBLE) THEN
      ! This should ideally be handled by the integrator (limits), 
      ! but returning the analytic limit is safest.
      integrand_result = current_mu - current_y
  
      RETURN
    ELSE
      CALL findRek(i, t, Rek)
      CALL findImk(i, t, Imk)
      
      integrand_result = DEXP( Rek ) * DSIN( Imk ) / t

    END IF
    
  END FUNCTION DFintegrand

END MODULE DFintegrand_MOD