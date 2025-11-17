MODULE Calcs_Solvers

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  PRIVATE
  
  ! Make the public subroutines available to other files using the module
  PUBLIC :: rtnewton, rtsafe, funcd_signature
  
  ! --- Interface Template ---
  ! This defines the signature that rtnewton and rtsafe rely on.
  INTERFACE
    SUBROUTINE funcd_signature(i, x, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(IN)    :: x
      INTEGER(C_INT), INTENT(IN)          :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: f, df
    END SUBROUTINE funcd_signature
  END INTERFACE

CONTAINS

  SUBROUTINE rtnewton(i, funcd, xstart, xacc, root)
    ! This function implements the Newton-Raphson method for finding a root
    ! of the function 'funcd' between bounds x1 and x2, starting at xstart.
    ! It includes a line-search safeguard to ensure x remains > 0.
    !
    ! funcd is a user-supplied subroutine which returns
    ! both the function value and the first derivative of the function.
    
    USE tweedie_params_mod, ONLY: Cverbose
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: xstart, xacc
    INTEGER(C_INT), INTENT(IN)        :: i
    REAL(KIND=C_DOUBLE)               :: root
    
    INTEGER, PARAMETER        :: MAXITS = 200
    INTEGER                   :: j
    REAL(KIND=C_DOUBLE)       :: dx, df, f
    REAL(KIND=C_DOUBLE)       :: x_current, x_iter_old, factor, x_new
    
    PROCEDURE(funcd_signature):: funcd
    
    
    ! Initialize starting guess
    x_current = xstart
    
    ! Check initial boundary condition
    IF (x_current .LE. 0.0_C_DOUBLE) THEN
      IF (Cverbose) CALL DBLEPR("ERROR (rtnewton): Initial guess is not positive:", -1, xstart, 1)
      root = HUGE(1.0_C_DOUBLE)
      RETURN
    END IF
    
    DO j = 1, MAXITS
      ! Save the old value, guaranteed to be > 0
      x_iter_old = x_current
      
      ! Calculate function value (f) and derivative (df) at the current safe point
      CALL funcd(i, x_current, f, df) 
      
      ! Check for convergence based on function value
      IF (ABS(f) .LT. xacc) THEN
        EXIT ! Root found (f is close to zero)
      END IF
      
      ! Check for near-zero derivative (Newton-Raphson failure)
      IF (ABS(df) .LT. 1.0E-12_C_DOUBLE) THEN ! Use a non-zero tolerance here
        IF (Cverbose) CALL DBLEPR("ERROR (rtnewton): Derivative near zero:", -1, df, 1)
        EXIT
      END IF
      
      ! Newton-Raphson Step: dx = f/f'
      dx = f / df
          
      ! Initialize step fraction and compute new value
      factor = 1.0_C_DOUBLE
      x_new = x_iter_old - dx ! Full step proposal
      
      ! SAFEGUARD: 
      ! Step halving (line search) to ensure x_new > 0
      DO WHILE ( (x_new .LE. 0.0_C_DOUBLE) .AND.   &
                 (factor .GT. 1.0E-6_C_DOUBLE) )
        ! If it violates the boundary, halve the step size
        factor = factor / 2.0_C_DOUBLE
        ! Re-calculate step based on the last safe value (x_iter_old)
        x_new = x_iter_old - factor * dx
      END DO  
      
      ! Check for catastrophic failure near zero
      IF (x_new .LE. 0.0_C_DOUBLE) THEN
          IF (Cverbose) CALL DBLEPR("ERROR (rtnewton): Step halving failed to find positive root:", -1, x_new, 1)
          ! If we can't step forward even a tiny bit, assume the root is effectively 0 or failed.
        x_current = x_iter_old
        EXIT
      END IF
    
      ! Update the current root estimate with the safe step
      x_current = x_new
      
      ! Check for step convergence (step size is small)
      ! Use the difference between the new and old value for robustness.
      IF (ABS(x_current - x_iter_old) .LT. xacc) EXIT
    END DO
    
    ! If we get here, NO CONVERGENCE after MAXITS
    IF (j .GE. MAXITS) THEN
      IF (Cverbose) CALL INTPR( "ERROR (rtnewton): failed to convergence after iterations:", -1, MAXITS, 1)
      IF (Cverbose) CALL DBLEPR("  Function value is:", -1, f, 1)
      IF (Cverbose) CALL DBLEPR("  at x:", -1, x_current, 1)
      IF (Cverbose) CALL DBLEPR("  and slope is:", -1, df, 1)
    END IF
    
    ! Assign the final root value
    root = x_current
  
  END SUBROUTINE rtnewton
  
  
  
  
  SUBROUTINE rtsafe(i, funcd, x1, x2, xacc, root)
    ! Adapted from NUMERICAL RECIPES Sect. 9.4
    ! Uses a combination of Newton-Raphson and bisection, find the root of a function bracketed
    ! between x1 and x2. The root is refined until its accuracy is known within plus/minus xacc. 
    !
    ! funcd is a user-supplied subroutine which returns
    ! both the function value and the first derivative of the function.
    
    USE tweedie_params_mod, ONLY: Cverbose
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
    
    IMPLICIT NONE
    
    INTEGER(C_INT), INTENT(IN)        :: i
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: x1, x2, xacc
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root
    
    INTEGER(C_INT)                    :: MAXIT, j
    REAL(KIND=C_DOUBLE)               :: df, dx, dxold, f, fh, fl, temp, xh, xl
    
    PARAMETER (MAXIT=100)   ! Maximum allowed number of iterations.
    
    PROCEDURE(funcd_signature) :: funcd
    
    call funcd(i, x1, fl, df)
    call funcd(i, x2, fh, df)
    
    IF ( (fl .GT. 0.0_C_DOUBLE) .AND. (fh .GT. 0.0_C_DOUBLE) &  
           .OR.                                                &
           (fl .LT. 0.0_C_DOUBLE) .AND. (fh .LT. 0.0_C_DOUBLE) ) THEN
      IF (Cverbose) CALL DBLEPR("ERROR (rtsafe): Root not bracketed: ", -1, x2, 1)
      STOP
    END IF
    
    IF (fl .EQ. 0.0_C_DOUBLE) THEN
      root = x1
      RETURN
    ELSE IF (fh .EQ. 0.0_C_DOUBLE) THEN
      root = x2
      RETURN
    ELSE IF (fl .LT. 0.0_C_DOUBLE) THEN 
    ! Orient the search so that f(xl) < 0.
      xl = x1
      xh = x2
    ELSE
      xh = x1
      xl = x2
    END IF
    
    root = 0.5_C_DOUBLE * (x1 + x2)   ! Initialize the guess for root,
    dxold = DABS(x2 - x1)             ! the “stepsize before last,”
    dx = dxold                        ! and the last step.
    CALL funcd(i, root, f, df)
    
    DO j = 1, MAXIT ! Loop over allowed iterations.
      IF ( ( (root - xh) * df - f) * ( (root - xl) * df - f) .GT. 0.0_C_DOUBLE    & 
             .OR.                                                                   &
             DABS(2.0_C_DOUBLE * f) .GT. DABS(dxold * df) ) THEN
        ! Bisect if Newton out of range,
        ! or not decreasing fast enough.
        dxold = dx
        dx = 0.5_C_DOUBLE *(xh - xl)
        root = xl + dx
        
        IF (xl .EQ. root) RETURN  ! Change in root is negligible.
      ELSE  ! Newton step acceptable; take it.
        dxold = dx
        dx = f/df
        temp = root
        root = root - dx
        IF (temp .EQ. root) RETURN
      END IF
      
      IF (DABS(dx) .LT. xacc) RETURN  ! Convergence criterion.
      CALL funcd(i, root, f, df)      ! The one new function evaluation per iteration.
      IF (f .LT. 0.0_C_DOUBLE) THEN   ! Maintain the bracket on the root.
        xl = root
      ELSE
        xh = root
      END IF
    END DO
    
    IF (j .GE. MAXIT) THEN
      IF (Cverbose) CALL INTPR( "ERROR (rtnewton): failed to convergence after iterations:", -1, MAXIT, 1)
      IF (Cverbose) CALL DBLEPR("  Function value is:", -1, f, 1)
      IF (Cverbose) CALL DBLEPR("  at x:", -1, root, 1)
      IF (Cverbose) CALL DBLEPR("  and slope is:", -1, df, 1)
    END IF
    
    RETURN
    
  END SUBROUTINE rtsafe

END MODULE Calcs_Solvers

