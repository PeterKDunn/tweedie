      FUNCTION rtnewton(funcd, x1, x2, xstart, xacc)
*     Using the Newton-Raphson method, find the root of a function known to lie in the interval
*     [x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
*     is a user-supplied subroutine that returns both the function value and the first derivative
*     of the function at the point x.
*
*     Used for find kmax and the zeros. All of these are positive,
*     and so we can found at zero on left.

      IMPLICIT NONE
      INTEGER JMAX, j
      DOUBLE PRECISION rtnewton, x1, x2, xacc, xstart
      DOUBLE PRECISION f, df, dx, rtnewtonTMP
      EXTERNAL funcd
      PARAMETER (JMAX=50) 
      
      write(*,*) "RTNEWTON:"
      rtnewton = xstart
    
      DO j = 1, JMAX
        CALL funcd(rtnewton, f, df)
        dx = f/df
        rtnewtonTMP = rtnewton - dx

*       Check that the new value does not suggest a negative value        
 50     IF (rtnewtonTMP .LT. 0.0d00) THEN
          write(*,*) "ADJUSTING for -ive value!"
          
          rtnewton = rtnewton / 2.0d00
          
          CALL funcd(rtnewton, f, df)
          dx = f/df
          rtnewtonTMP = rtnewton - dx

          GOTO 50
        ENDIF
        rtnewton = rtnewton
*      write(*,*) "       j:", j
*      write(*,*) "rtnewton:", rtnewton
*      write(*,*) "       f:", f
*      write(*,*) "      df:", df
*      write(*,*) "  error:", DABS(dx)
*      write(*,*) "   xacc:", xacc

        IF ( (x1 - rtnewton) * (rtnewton - x2) .LT. 0.0d00) THEN
*          write(*,*) "rtnewton:", rtnewton
*          write(*,*) "       f:", f
*          write(*,*) "      df:", df
*          write(*,*) "rtnewton jumped out of brackets"
*          STOP
           rtnewton = DABS(rtnewton)
*      write(*,*) "   -> rtnewton:", rtnewton
        ENDIF
        IF ( DABS(dx) .LT. xacc) RETURN 
*       Convergence!
      ENDDO
*      write(*,*) "rtnewton: exceeded maximum iterations", f
*      write(*,*) "       j:", j
*      write(*,*) "rtnewton:", rtnewton
*      write(*,*) "       f:", f
*      write(*,*) "      df:", df
      
      RETURN
      END