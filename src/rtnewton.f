      FUNCTION rtnewton(funcd, x1, x2, xstart, xacc)
*     Using the Newton-Raphson method, find the root of a function known to lie in the interval
*     [x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
*     is a user-supplied subroutine that returns both the function value and the first derivative
*     of the function at the point x.
        
      IMPLICIT NONE
      INTEGER JMAX, j
      DOUBLE PRECISION rtnewton, x1, x2, xacc, xstart
      DOUBLE PRECISION f, df, dx
      EXTERNAL funcd
      PARAMETER (JMAX=50) 
      
      rtnewton = xstart
        
      DO j = 1, JMAX
        CALL funcd(rtnewton, f, df)
        dx = f/df
        rtnewton = rtnewton - dx

        IF ( (x1 - rtnewton) * (rtnewton - x2) .LT. 0.0d00) THEN
          write(*,*) "rtnewton:", rtnewton
          write(*,*) "       f:", f
          write(*,*) "      df:", df
          write(*,*) "rtnewton jumped out of brackets"
*          STOP
        ENDIF
        IF ( DABS(dx) .LT. xacc) RETURN 
*       Convergence!
      ENDDO
      write(*,*) "rtnewton: exceeded maximum iterations"
      write(*,*) "rtnewton:", rtnewton
      write(*,*) "       f:", f
      write(*,*) "      df:", df
      STOP
      
      RETURN
      END