SUBROUTINE rtnewton(i, funcd, xstart, xacc, root)
  ! This function implements the Newton-Raphson method for finding a root
  ! of the function 'funcd' between bounds x1 and x2, starting at xstart.
  ! It includes a line-search safeguard to ensure x remains > 0.
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

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

  PROCEDURE(funcd_signature):: funcd
  

  REAL(KIND=C_DOUBLE), INTENT(IN)  :: xstart, xacc
  INTEGER(C_INT), INTENT(IN)        :: i
  
  ! Output (Function result)
  REAL(KIND=C_DOUBLE)              :: root
  
  ! --- Local Variables
  INTEGER, PARAMETER  :: MAXITS = 5000
  INTEGER             :: j
  
  ! x_current holds the final accepted value, x_iter_old is the last safe value.
  ! x_new holds the proposed step value.
  REAL(KIND=C_DOUBLE)      :: dx, df, f
  REAL(KIND=C_DOUBLE)      :: x_current, x_iter_old, factor, x_new
  
  ! Initialize starting guess
  x_current = xstart

  ! Check initial boundary condition
  IF (x_current <= 0.0_C_DOUBLE) THEN
      WRITE(*,*) "RTNEWTON ERROR: Initial guess (xstart) is not positive."
      root = HUGE(1.0_C_DOUBLE)
      RETURN
  END IF
  
  DO j = 1, MAXITS

    ! Save the old value, which is guaranteed to be > 0
    x_iter_old = x_current
    
    ! 1. Calculate function value (f) and derivative (df) at the current safe point
    CALL funcd(i, x_current, f, df) 

    ! 2. Check for convergence based on function value
    IF (ABS(f) .LT. xacc) THEN
        EXIT ! Root found (f is close to zero)
    END IF
    
    ! 3. Check for near-zero derivative (Newton-Raphson failure)
    IF (ABS(df) .LT. 1.0E-12_C_DOUBLE) THEN ! Use a non-zero tolerance here
      WRITE(*,*) "RTNEWTON ERROR: Derivative near zero."
      EXIT
    END IF
    
    ! 4. Newton-Raphson Step: dx = f/f'
    dx = f / df
    
    ! 5. Initialize step fraction and proposed new value
    factor = 1.0_C_DOUBLE
    x_new = x_iter_old - dx ! Full step proposal
    
    ! 6. Safeguard Check: Step Halving (Line Search) to ensure x_new > 0
    DO WHILE (x_new <= 0.0_C_DOUBLE .AND. factor > 1.0E-6_C_DOUBLE)
        ! If it violates the boundary, halve the step size
        factor = factor / 2.0_C_DOUBLE
        ! Re-calculate step based on the last safe value (x_iter_old)
        x_new = x_iter_old - factor * dx
    END DO  
    
    ! 7. Check for catastrophic failure near zero
    IF (x_new <= 0.0_C_DOUBLE) THEN
        WRITE(*,*) "RTNEWTON ERROR: Step halving failed to find positive x."
        ! If we can't step forward even a tiny bit, assume the root is effectively 0 or failed.
        x_current = x_iter_old
        EXIT
    END IF

    ! 8. Update the current root estimate with the safe step
    x_current = x_new

    ! 9. Check for step convergence (step size is small)
    ! We use the difference between the new and old value for robustness.
    IF (ABS(x_current - x_iter_old) < xacc) EXIT
  END DO
  
  ! If we get here, NO CONVERGENCE after MAXITS
  IF (j .GE. MAXITS) THEN
    WRITE(*,*) "CONVERGENCE NOT OBTINED IN RTNEWTON after", MAXITS, "iterations"
    WRITE(*,*) " - Have f, df = ", f, df
    
  END IF

  ! Assign the final root value
  root = x_current

END SUBROUTINE rtnewton
