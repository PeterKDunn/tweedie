SUBROUTINE rtsafe(i, funcd, x1, x2, xstart, xacc, root)
  USE tweedie_params_mod

  IMPLICIT NONE
  
  INTERFACE
      SUBROUTINE funcd_signature(i, x, f, df)
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

        IMPLICIT NONE
        REAL(KIND=C_DOUBLE), INTENT(IN)    :: x
        INTEGER(C_INT), INTENT(IN)         :: i
        REAL(KIND=C_DOUBLE), INTENT(OUT)   :: f, df

      END SUBROUTINE funcd_signature
  END INTERFACE
  
  PROCEDURE(funcd_signature) :: funcd
 

  ! Inputs
  REAL(KIND=C_DOUBLE), INTENT(IN)  :: x1, x2, xstart, xacc
  INTEGER(C_INT), INTENT(IN)       :: i

  ! Output (Function result)
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: root
  
  INTEGER, PARAMETER   :: JMAX = 50
  INTEGER              :: j
  REAL(KIND=8)         :: f, df, dx, rootTMP
  REAL(KIND=8)         :: rtsafeTMP, fh, fl, xh, xl
  
  ! Initialize root and boundaries
  xl = x1
  xh = x2
  rootTMP = xstart
  
  ! Get function values at boundaries
  CALL funcd(i, xl, fl, df)
  CALL funcd(i, xh, fh, df)
  
  ! --- Main Loop ---
  
  DO j = 1, JMAX
    CALL funcd(i, rootTMP, f, df)
    
    ! Check for non-finite derivative (NaN/Inf) OR very small derivative
    IF ( (df .NE. df) .OR. (DABS(df) < 1.0E-12_C_DOUBLE) ) THEN
        ! Fallback: Use bisection if Newton's method is unstable
        rtsafeTMP = (xl + xh) / 2.0_C_DOUBLE
        dx = DABS(rtsafeTMP - rootTMP)
        rootTMP = rtsafeTMP
    ELSE
        ! Use Newton-Raphson step
        dx = f / df
        rtsafeTMP = rootTMP - dx
        
        ! Check if Newton step is valid (remains within brackets)
        IF (((rootTMP - rtsafeTMP) * f) .LT. 0.0_C_DOUBLE) THEN
            rootTMP = rtsafeTMP
        ELSE
            ! Fallback to bisection if step is invalid
            rtsafeTMP = (xl + xh) / 2.0_C_DOUBLE
            dx = DABS(rtsafeTMP - rootTMP)
            rootTMP = rtsafeTMP
        END IF
    END IF
    
    ! Check convergence (using the step size from the method used)
    IF (DABS(dx) < xacc) RETURN 
    
    ! Update brackets (bisection logic)
    IF (f .LT. 0.0_C_DOUBLE) THEN
        xl = rootTMP
    ELSE
        xh = rootTMP
    END IF
    
  END DO
  
  root = rootTMP
  
END SUBROUTINE rtsafe
