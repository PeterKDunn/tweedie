! Eliminates the F77 COMMON block and updates calls to findImk and findRek
! to prevent segmentation faults.
!
FUNCTION DFintegrand(t) RESULT(integrand_result) BIND(C, NAME='DFintegrand')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
 ! NOTE: Added INTERFACE blocks to declare C-binding for helper routines
  INTERFACE
    SUBROUTINE findImk(i, t, Imk) BIND(C, NAME='findImk')
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      
      INTEGER(C_INT), INTENT(IN)          :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Imk
    END SUBROUTINE findImk
    
    SUBROUTINE findRek(i, t, Rek) BIND(C, NAME='findRek')
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)          :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Rek
    END SUBROUTINE findRek
  END INTERFACE
    
  ! Arguments (Inputs)
  REAL(KIND=C_DOUBLE), INTENT(IN) :: t                ! The internal variable for integration
  REAL(KIND=C_DOUBLE)             :: integrand_result ! The result of the function

  REAL(KIND=8)          :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE)   :: Imk, Rek


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(current_i)
  current_mu   = Cmu(current_i)
  current_phi  = Cphi(current_i)
  
  ! The call to findImk and findRek must pass the parameters explicitly.
  
  ! Check for when t = 0 (t is very close to zero)
  IF (ABS(t) < 1.0d-14) THEN
    ! This should ideally be handled by the integrator (limits), 
    ! but returning the analytic limit is safest.
    integrand_result = current_mu - current_y

    RETURN
  ELSE
    ! NOTE: Calls to findImk and findRek MUST include the extra parameters
    CALL findImk(current_i, t, Imk)
    CALL findRek(current_i, t, Rek) ! Note: Cy passed to findRek
    
    integrand_result = EXP( Rek ) * SIN( Imk ) / t

  END IF
  
END FUNCTION DFintegrand
