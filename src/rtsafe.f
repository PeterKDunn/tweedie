      FUNCTION rtsafe(funcd, x1, x2, xstart, xacc)
*     Newton-Raphson root finder with bounds safeguard
*     Ensures root stays within [x1, x2] and never goes negative
      IMPLICIT NONE
      INTEGER, PARAMETER :: JMAX = 50
      INTEGER :: j
      DOUBLE PRECISION :: rtsafe, x1, x2, xstart, xacc
      DOUBLE PRECISION :: f, df, dx, rtsafeTMP
      EXTERNAL funcd

      ! Initialize root
      rtsafe = xstart

*      write(*,*) "RTSAFE: Bounded between ", x1, x2
      DO j = 1, JMAX
         CALL funcd(rtsafe, f, df)
         dx = f / df
         rtsafeTMP = rtsafe - dx

*        Keep root within [x1, x2]
         IF (rtsafeTMP < x1) THEN
            rtsafeTMP = 0.5D0 * (rtsafe + x1)
         ELSE IF (rtsafeTMP > x2) THEN
            rtsafeTMP = 0.5D0 * (rtsafe + x2)
         ENDIF

         rtsafe = rtsafeTMP
*      write(*,*) "RTSAFE: x, fn value: ", rtsafe, f

*        Convergence check
         IF (DABS(dx) < xacc) RETURN
      ENDDO

      ! Optional: warn if maximum iterations exceeded
      WRITE(*,*) 'rtsafe: WARNING - maximum iterations reached'

      RETURN
      END
