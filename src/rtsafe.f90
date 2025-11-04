SUBROUTINE rtsafe(i, funcd, x1, x2, xacc, root)
 ! Adapted from NUMERICAL RECIPES Sect. 9.4
!  FUNCTION rtsafe(funcd, x1, x2, xacc)

  USE tweedie_params_mod, ONLY: Cverbose
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTERFACE
      SUBROUTINE funcd_signature(i, x, f, df)
      ! Template for the function for which roots are sought.
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

        IMPLICIT NONE
        INTEGER(C_INT), INTENT(IN)          :: i
        REAL(KIND=C_DOUBLE), INTENT(IN)     :: x
        REAL(KIND=C_DOUBLE), INTENT(OUT)    :: f, df
      END SUBROUTINE funcd_signature
  END INTERFACE

  PROCEDURE(funcd_signature) :: funcd

  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: x1, x2, xacc
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root

  INTEGER(C_INT)                    :: MAXIT, j
  REAL(KIND=C_DOUBLE)               :: df, dx, dxold, f, fh, fl, temp, xh, xl
  
  PARAMETER (MAXIT=100)   ! Maximum allowed number of iterations.
   ! Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
  ! between x1 and x2. The root, returned as the function value rtsafe, will be refined until
  ! its accuracy is known within ±xacc. funcd is a user-supplied subroutine which returns
  ! both the function value and the first derivative of the function.



  call funcd(i, x1, fl, df)
  call funcd(i, x2, fh, df)
  
  IF ( (fl .GT. 0.0_C_DOUBLE) .AND. (fh .GT. 0.0_C_DOUBLE) &  
       .OR.                                                &
       (fl .LT. 0.0_C_DOUBLE) .AND. (fh .LT. 0.0_C_DOUBLE) ) THEN
    WRITE(*,*) "ERROR (rtsafe): Root must be bracketed in rtsafe"
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

  root = 0.5_C_DOUBLE * (x1 + x2) ! Initialize the guess for root,
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
    ELSE  ! Newton step acceptable. Take it.
      dxold = dx
      dx = f/df
      temp = root
      root = root - dx
      IF (temp .EQ. root) RETURN
    END IF
  
    IF (DABS(dx) .LT. xacc) RETURN ! Convergence criterion.
    CALL funcd(i, root, f, df) ! The one new function evaluation per iteration.
    IF (f .LT. 0.0_C_DOUBLE) THEN ! Maintain the bracket on the root.
      xl = root
    ELSE
      xh = root
    END IF
  END DO

  IF (Cverbose) WRITE (*,*) "Perhaps non-convergence (rtsafe): exceeding maximum iterations"
  RETURN

END SUBROUTINE rtsafe
