MODULE Calcs_Real

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE R_interfaces

  IMPLICIT NONE
  
  PUBLIC :: evaluateRek, evaluateRekd
  
CONTAINS

  SUBROUTINE evaluateRek(t, Rek, errorHere)
    ! Find the value of Re k(t)
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Rek
    LOGICAL(C_BOOL), INTENT(OUT)        :: errorHere    ! TRUE if any computational problem founds

    REAL(KIND=C_DOUBLE) :: omega, pindex, front, alpha, tanArg, zRe, zIm
  
    
    errorHere = .FALSE.
    Rek       = 0.0E0_C_DOUBLE
    
    pindex = (2.0E0_C_DOUBLE - Cp)
    front  = current_mu ** pindex  / ( current_phi * pindex)
    tanArg = (1.0E0_C_DOUBLE - Cp) * t * current_phi / (current_mu ** (1.0E0_C_DOUBLE - Cp) )
    omega  = DATAN( tanArg )
    
    ! Safety check
    IF ((omega .GT. 0.0E0_C_DOUBLE ) .OR. (omega .LT. (-PI/2.0E0_C_DOUBLE))) THEN
      ! Error!
      errorHere = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: Rek: omega out of bounds =", -1, omega, 1)
      
      RETURN
    END IF
    
    alpha = (2.0E0_C_DOUBLE - Cp)/(1.0E0_C_DOUBLE - Cp)
    ! Re k(t) = front * ( Re[(1 + i tanArg)^alpha] - 1 )   (stably: see zPow)
    CALL zPow(tanArg, alpha, zRe, zIm)
    Rek = front * ( zRe - 1.0E0_C_DOUBLE )

  END SUBROUTINE evaluateRek



  SUBROUTINE evaluateRekd(t, Rekd, errorHere)
    ! Find the value of Re k'(t)
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Rekd
    LOGICAL(C_BOOL), INTENT(OUT)      :: errorHere
    REAL(KIND=C_DOUBLE)               :: omega, pindex, tanArg, zRe, zIm

    ! Initialise
    errorHere = .FALSE.
    Rekd      = 0.0E0_C_DOUBLE
    
    pindex = 1.0E0_C_DOUBLE / (1.0E0_C_DOUBLE - Cp)
    tanArg = ( (1.0E0_C_DOUBLE - Cp) * t * current_phi) / &
             (current_mu ** (1.0E0_C_DOUBLE - Cp) )
    omega  = DATAN( tanArg )
    
    ! Safety check
    IF ((omega .GT. 0.0E0_C_DOUBLE ) .OR. (omega .LT. (-PI/2.0E0_C_DOUBLE))) THEN
      ! Error!
      errorHere = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: Rekd: omega out of bounds =", -1, omega, 1)
      RETURN
    END IF

    ! Re k'(t) = mu * Im[(1 + i tanArg)^pindex]   (stably: see zPow)
    CALL zPow(tanArg, pindex, zRe, zIm)
    Rekd = current_mu * zIm
    
  END SUBROUTINE evaluateRekd
  
  
  
  
    
  SUBROUTINE evaluateLambda(lambda)
    ! Find lambda, such that P(Y = 0) = exp( -lambda ) when 1 < p < 2 
    
    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: lambda 
    

    lambda = 0.0E0_C_DOUBLE
    IF (CpSmall) THEN
      ! The calculation for lambda (used in P(Y=0) = exp(-lambda))
      lambda = (current_mu ** (2.0E0_C_DOUBLE - Cp) ) / &
               (current_phi * (2.0E0_C_DOUBLE - Cp) )
      ! NOTE: No negative sign in front
    END IF
    
  END SUBROUTINE evaluateLambda

END MODULE Calcs_Real
