SUBROUTINE rtnewton(i, funcd, x1, x2, xstart, xacc, root)
  ! This function implements the Newton-Raphson method for finding a root
  ! of the function 'funcd' between bounds x1 and x2, starting at xstart.
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
  

  REAL(KIND=C_DOUBLE), INTENT(IN)  :: x1, x2, xstart, xacc
  INTEGER(C_INT), INTENT(IN)       :: i
  
  ! Output (Function result)
  REAL(KIND=C_DOUBLE)              :: root
  
  ! --- Local Variables
  INTEGER, PARAMETER  :: MAXITS = 5000
  INTEGER             :: j
  
  ! x_current holds the current iteration value, dx is the change
  REAL(KIND=C_DOUBLE)        :: dx, df, f, x_current, x_iter_old
  

  x_current = xstart ! Start at the initial guess
  
  DO j = 1, MAXITS
!  write(*,*) " IN rtnewton; j=",j
!  write(*,*), "   with x:", x_current
    ! Save the old value in case we need to revert
    x_iter_old = x_current
    
    ! Calculate function value (f) and derivative (df)
    CALL funcd(i, x_current, f, df) 

!  write(*,*), "   with f, df:",f, df 
    ! Check for convergence before a step
    IF (ABS(f) < xacc) THEN
        EXIT ! Root found (f is close to zero)
    END IF
    
    IF (ABS(df) < xacc) THEN
      ! Derivative is near zero: Newton-Raphson fails. 
      ! Revert to last safe value and exit to force error/fallback (if required)
      x_current = x_iter_old
      EXIT
    END IF
    
    ! Newton-Raphson Step: dx = f/f'
    dx = f / df
    x_current = x_current - dx
    
    ! Check for bounds violation
    IF ( (x_current < x1) .OR. (x_current > x2) ) THEN
        ! Root is outside bounds, revert to last safe value and exit
        x_current = x_iter_old
        EXIT
    END IF
    
    ! Check for step convergence (step size is small)
    IF (ABS(dx) < xacc) EXIT
  END DO
  
  ! If we get here, NO CONVERGENCE
  ! Assign the final root value to the function result
  write(*,*) "CONVERGENCE NOT OBTINED IN RTNEWTON; t, f, df:", x_current, f, df
  root = x_current

END SUBROUTINE rtnewton
