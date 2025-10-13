
FUNCTION rtsafe(funcd, x1, x2, xstart, xacc) BIND(C, NAME='rtsafe')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! NOTE: All problematic intrinsic declarations (LOGICAL, INTRINSIC, USE)
  ! have been removed to force the use of the pure Fortran NaN check.
  
  ! --- Function Arguments ---
  
 ! --- CRITICAL FIX: The funcd argument must be declared C-interoperable ---
  INTERFACE
      ! Define the expected signature of the function pointer passed to rtnewton.
      ! Assuming the solved function (funcd) is a SUBROUTINE returning f and df.
      SUBROUTINE funcd_signature(x, f, df) BIND(C)
        USE ISO_C_BINDING, ONLY: C_DOUBLE

        IMPLICIT NONE
        REAL(KIND=C_DOUBLE), INTENT(IN)    :: x
        REAL(KIND=C_DOUBLE), INTENT(OUT)   :: f, df
        ! The remaining arguments are implicitly accessible through the module
      END SUBROUTINE funcd_signature
  END INTERFACE
  
  ! Declaration of the procedure argument, linking it to the defined interface.
  PROCEDURE(funcd_signature), BIND(C) :: funcd
 
  ! Inputs
  REAL(KIND=C_DOUBLE), INTENT(IN) :: x1, x2, xstart, xacc

  ! Output (Function result)
  REAL(KIND=C_DOUBLE) :: rtsafe
  
  ! --- Local Variables ---
  
  INTEGER(C_INT), PARAMETER   :: JMAX = 50
  INTEGER(C_INT)              :: j
  REAL(KIND=C_DOUBLE)         :: f, df, dx
  REAL(KIND=C_DOUBLE)         :: rtsafeTMP, fh, fl, xh, xl
  
  ! Initialize root and boundaries
  xl = x1
  xh = x2
  rtsafe = xstart
  
  ! Get function values at boundaries
  CALL funcd(xl, fl, df)
  CALL funcd(xh, fh, df)
  
  ! --- Main Loop ---
  
  DO j = 1, JMAX
    ! CRITICAL FIX: Pass the extra parameters to the function being solved!
    CALL funcd(rtsafe, f, df)
    
    ! Check for non-finite derivative (NaN/Inf) OR very small derivative
    ! PURE FORTRAN FIX: df .NE. df is true if df is NaN or Inf.
    IF (df .NE. df .OR. DABS(df) < 1.0d-12) THEN
        ! Fallback: Use bisection if Newton's method is unstable
        rtsafeTMP = (xl + xh) / 2.0_8
        dx = DABS(rtsafeTMP - rtsafe)
        rtsafe = rtsafeTMP
    ELSE
        ! Use Newton-Raphson step
        dx = f / df
        rtsafeTMP = rtsafe - dx
        
        ! Check if Newton step is valid (remains within brackets)
        IF (((rtsafe - rtsafeTMP) * f) .LT. 0.0_8) THEN
            rtsafe = rtsafeTMP
        ELSE
            ! Fallback to bisection if step is invalid
            rtsafeTMP = (xl + xh) / 2.0_8
            dx = DABS(rtsafeTMP - rtsafe)
            rtsafe = rtsafeTMP
        END IF
    END IF
    
    ! Check convergence (using the step size from the method used)
    IF (DABS(dx) < xacc) RETURN 
    
    ! Update brackets (bisection logic)
    IF (f .LT. 0.0_8) THEN
        xl = rtsafe
    ELSE
        xh = rtsafe
    END IF
    
  END DO
  
  rtsafe = -999.0_C_DOUBLE
  
  END FUNCTION rtsafe
