SUBROUTINE rtsafe(i, funcd, xstart, x1, x2, xacc, root)
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
  
  PROCEDURE(funcd_signature) :: funcd
 

  INTEGER(C_INT), INTENT(IN)       :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)  :: x1, x2, xstart, xacc
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: root
  
  INTEGER, PARAMETER    :: JMAX = 500
  INTEGER               :: j
  REAL(KIND=C_DOUBLE)   :: f, df, dx, rootTMP
  REAL(KIND=C_DOUBLE)   :: rtsafeTMP, fh, fl, xh, xl
  
  ! Initialize root and boundaries
  xl = x1
  xh = x2
  rootTMP = xstart
  
  ! Get function values at boundaries
  CALL funcd(i, xl, fl, df)
  CALL funcd(i, xh, fh, df)
  
  ! --- Main Loop ---
  
  IF (fl * fh > 0.0_C_DOUBLE) THEN
      ! Error: The bounds do not bracket a root. 
      ! This should halt execution or return an error flag.
      WRITE(*,*) "RTSAFE ERROR: ROOT NOT BRACKETED BY X1 AND X2"
      WRITE(*,*) "RTSAFE ERROR: xl, xh", xl, xh
      WRITE(*,*) "RTSAFE ERROR: fl, fh", fl, fh
      root = HUGE(1.0_C_DOUBLE) ! Return a clear error value
      RETURN 
  ELSE IF (fl > 0.0_C_DOUBLE) THEN
      ! If fl is positive, swap xl and xh (and fl and fh) to ensure fl < 0
      rootTMP = xl
      xl = xh
      xh = rootTMP
      rootTMP = fl
      fl = fh
      fh = rootTMP
  END IF

! At this point, we guarantee: fl < 0 and fh > 0.
WRITE(*,*) "STARTING:: x", xstart


  DO j = 1, JMAX    
      ! 1. Evaluate function and derivative at current estimate (rootTMP)
      CALL funcd(i, rootTMP, f, df)
      
      ! Check for convergence based on function value (recommended)
      IF (DABS(f) < xacc) EXIT 

      ! 2. Compute proposed Newton step
      ! This is ONLY the proposed value, not yet the accepted value
      dx = f / df
      rtsafeTMP = rootTMP - dx
      
      ! -----------------------------------------------------------
      ! 3. SAFEGUARD CHECK: DECIDE BETWEEN NEWTON AND BISECTION
      ! -----------------------------------------------------------
      
      ! Check 1: Is the derivative flat (ABS(df) is NOT too small)? AND
      ! Check 2: Does the Newton step stay within the current, shrinking bracket [xl, xh]?
! -----------------------------------------------------------
      ! 3. SAFEGUARD CHECK: DECIDE BETWEEN NEWTON AND BISECTION
      ! -----------------------------------------------------------
      
      ! Use Newton Step if it lands within the CURRENT safe bracket AND the derivative is not flat.
      IF (rtsafeTMP > xl .AND. rtsafeTMP < xh .AND. DABS(df) > 1.0E-12_C_DOUBLE) THEN
          ! Newton Step is SAFE and accepted.
      ELSE
          ! Newton Step is rejected (overshot the bracket OR derivative is flat).
          ! Fallback to Bisection Step
          rtsafeTMP = (xl + xh) / 2.0_C_DOUBLE
      END IF      
      ! -----------------------------------------------------------
      ! 4. Update the root and check step size convergence
      ! -----------------------------------------------------------
      
      ! The step size is the difference between the new accepted value and the old value
      dx = rtsafeTMP - rootTMP 
      rootTMP = rtsafeTMP
      
      ! Check convergence based on step size
      IF (DABS(dx) < xacc) EXIT
      
      ! -----------------------------------------------------------
      ! 5. SHRINK BRACKET: Update bounds based on the sign of f
      ! -----------------------------------------------------------
      
      ! Re-evaluate f at the newly accepted rootTMP to update the brackets correctly
      ! (We must re-evaluate because the initial funcd call was at the old rootTMP)
      CALL funcd(i, rootTMP, f, df) 
      
      ! Shrink the bracket [xl, xh] based on the sign of f at rootTMP
      ! This requires that the function was initially bracketed with f(xl)*f(xh) < 0
      IF (f < 0.0_C_DOUBLE) THEN
          xl = rootTMP
      ELSE
          xh = rootTMP
      END IF
      
      ! Check for bracket width convergence (recommended)
      IF (DABS(xh - xl) < xacc) EXIT
WRITE(*,*) "j, x, f, df", j, rootTMP, f, df

  END DO
root = rootTMP
STOP
  
END SUBROUTINE rtsafe
