      FUNCTION rtnewton(funcd, x1, x2, xstart, xacc)
*     Using the Newton-Raphson method, find the root of a function known to lie in the interval
*     [x1, x2]. The root rtnewt will be refined until its accuracy is known within xacc. funcd
*     is a user-supplied subroutine that returns both the function value and the first derivative
*     of the function at the point x.
        
      IMPLICIT NONE
      INTEGER JMAX, j
      DOUBLE PRECISION rtnewton, x1, x2, xacc, xstart
      DOUBLE PRECISION f, df, dx
      EXTERNAL funcd



      JMAX = 50
      rtnewton = xstart
        
      DO j = 1, JMAX

        CALL funcd(rtnewton, f, df)
        dx = f / df
        rtnewton = rtnewton - dx
        
        IF ( (x1 - rtnewton) * (rtnewton - x2) .LT. 0.0d00) THEN
          IF ( rtnewton .LT. 0.0d00 ) THEN
            write(*,*) ">>> HACK: absolute value-ing"
            rtnewton = DABS(rtnewton)
          ENDIF
*          write(*,*) "rtnewton:", rtnewton
*          write(*,*) "       f:", f
*          write(*,*) "      df:", df
*          write(*,*) "rtnewton jumped out of brackets"
*          write(*,*) "These were", x1, x2
*          write(*,*) "Latest guess", rtnewton
*         ONLY use Newtons when we can be sure it converges,
*         which means we don;t have to be too concerned about this.
          x2 = x2 * 10.0d00
          x1 = x1 / 10.0d00
*          STOP
        ENDIF
        
        IF ( DABS(dx) .LT. xacc) THEN
*         Convergence!
*          write(*,*) "rtnewton: t =", rtnewton
*          write(*,*) "          f =", f
          RETURN 
        ENDIF
      ENDDO
      write(*,*) "rtnewton: exceeded maximum iterations"

      RETURN
      END