SUBROUTINE rtsafe(i, funcd, xstart, x1, x2, xacc, root)
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  IMPLICIT NONE

  INTERFACE
      SUBROUTINE funcd_signature(i, x, f, df)
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
        IMPLICIT NONE
        INTEGER(C_INT), INTENT(IN)         :: i
        REAL(KIND=C_DOUBLE), INTENT(IN)    :: x
        REAL(KIND=C_DOUBLE), INTENT(OUT)   :: f, df
      END SUBROUTINE funcd_signature
  END INTERFACE

  PROCEDURE(funcd_signature) :: funcd

  INTEGER(C_INT), INTENT(IN)       :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)  :: x1, x2, xstart, xacc
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: root

  INTEGER, PARAMETER    :: JMAX = 500
  INTEGER               :: j
  REAL(KIND=C_DOUBLE)   :: f, df, fl, fh
  REAL(KIND=C_DOUBLE)   :: xl, xh, dx, rtsafeTMP
  REAL(KIND=C_DOUBLE)   :: fold, testtol
  REAL(KIND=C_DOUBLE)   :: halfbracket, maxstep

  ! Initialize
  xl = x1
  xh = x2
  rtsafeTMP = xstart

  ! If xstart is outside the bracket, use midpoint
  IF ( (rtsafeTMP .LT. MIN(x1, x2) ) .OR.  &
       (rtsafeTMP .GT. MAX(x1, x2)) ) THEN
      rtsafeTMP = 0.5_C_DOUBLE * (x1 + x2)
  END IF

  ! Evaluate function at endpoints
  CALL funcd(i, xl, fl, df)
  CALL funcd(i, xh, fh, df)

  ! Check bracketing
  IF (fl * fh .GT. 0.0_C_DOUBLE) THEN
      WRITE(*,*) "RTSAFE ERROR: ROOT NOT BRACKETED BY X1 AND X2"
      WRITE(*,*) "  xl, xh = ", xl, xh
      WRITE(*,*) "  fl, fh = ", fl, fh
      root = HUGE(1.0_C_DOUBLE)
      RETURN
  END IF

  ! Ensure fl < 0 < fh (swap if necessary)
  IF (fl .GT. 0.0_C_DOUBLE) THEN
      ! swap xl <-> xh and fl <-> fh using temporaries (no separate subroutine)
      dx = xl
      xl = xh
      xh = dx
      dx = fl
      fl = fh
      fh = dx
  END IF

  ! Main loop
  DO j = 1, JMAX
      ! Evaluate function and derivative at current estimate
      CALL funcd(i, rtsafeTMP, f, df)

      ! Convergence on function value
      IF (ABS(f) .LT. xacc) THEN
          root = rtsafeTMP
          RETURN
      END IF

      ! Compute half-bracket width & limit for step
      halfbracket = 0.5_C_DOUBLE * (xh - xl)
      maxstep = MAX(ABS(rtsafeTMP), 1.0_C_DOUBLE) * xacc  ! relative related scale

      ! Decide whether to use Newton step or bisection:
      ! - if df is too small (relative), or Newton would step outside bracket,
      !   or Newton step would be unreasonably large -> use bisection (midpoint)
      IF (ABS(df) .LT. (1.0E-16_C_DOUBLE * (1.0_C_DOUBLE + ABS(f))) ) THEN
          ! derivative effectively zero -> midpoint
          rtsafeTMP = 0.5_C_DOUBLE * (xl + xh)
      ELSE
          ! proposed Newton step
          dx = f / df
          ! Limit Newton step so it does not exceed half the bracket
          IF (ABS(dx) .GT. halfbracket) THEN
              dx = SIGN(halfbracket, dx)
          END IF
          rtsafeTMP = rtsafeTMP - dx
          ! if Newton step jumps outside bracket, fall back to midpoint
          IF ( (rtsafeTMP .LE. xl) .OR. &
               (rtsafeTMP .GE. xh) ) THEN
              rtsafeTMP = 0.5_C_DOUBLE * (xl + xh)
          END IF
      END IF

      ! Evaluate function at accepted point once
      CALL funcd(i, rtsafeTMP, f, df)

      ! Shrink the bracket using the new f
      IF (f .LT. 0.0_C_DOUBLE) THEN
          xl = rtsafeTMP
          fl = f
      ELSE
          xh = rtsafeTMP
          fh = f
      END IF

      ! Convergence based on interval width (relative) or absolute small step
      testtol = xacc * (1.0_C_DOUBLE + ABS(rtsafeTMP))
      IF ((xh - xl) .LE. testtol) THEN
          root = rtsafeTMP
          RETURN
      END IF

  END DO

  ! If we exit loop without returning, max iterations reached -> return last estimate
  WRITE(*,*) "RTSAFE WARNING: maximum iterations reached without full convergence."
  root = rtsafeTMP

END SUBROUTINE rtsafe
