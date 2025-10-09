      FUNCTION rtnewton(funcd, x1, x2, xstart, xacc)
*     Drop-in replacement for rtnewtonOLD
*     Ensures the value never goes negative
      IMPLICIT NONE
      INTEGER JMAX, j
      DOUBLE PRECISION rtnewton, x1, x2, xstart, xacc
      DOUBLE PRECISION f, df, dx, rtnewtonTMP
      EXTERNAL funcd
      PARAMETER (JMAX=50)

      rtnewton = xstart

      DO j = 1, JMAX
        CALL funcd(rtnewton, f, df)
        dx = f / df
        rtnewtonTMP = rtnewton - dx
*        write(*,*) "RTNEWTON"
*        write(*,*) rtnewton, f, df
*       Safeguard: do not let root go negative
 50     IF (rtnewtonTMP .LT. 0.0D0) THEN
          rtnewton = rtnewton / 2.0D0
          CALL funcd(rtnewton, f, df)
          dx = f / df
          rtnewtonTMP = rtnewton - dx
          GOTO 50
        ENDIF

        rtnewton = rtnewtonTMP

*       Convergence check
        IF (DABS(dx) .LT. xacc) RETURN
      ENDDO

      RETURN
      END
