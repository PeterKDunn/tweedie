      FUNCTION rtnewtonOLD(funcd, x1, x2, xstart, xacc)
*     Using the Newton-Raphson method, find the root of a function known to lie in the interval
*     [x1, x2]. The root rtnewt will be refined until its accuracy is known within plus/minus xacc. funcd
*     is a user-supplied subroutine that returns both the function value and the first derivative
*     of the function at the point x.
*
*     Used for find kmax and the zeros. All of these are positive,
*     and so we can found at zero on left.

      IMPLICIT NONE
      INTEGER JMAX, j
      DOUBLE PRECISION rtnewtonOLD, x1, x2, xacc, xstart
      DOUBLE PRECISION f, df, dx, rtnewtonOLDTMP
      EXTERNAL funcd
      PARAMETER (JMAX=50) 
      
      rtnewtonOLD = xstart
    
      DO j = 1, JMAX
        CALL funcd(rtnewtonOLD, f, df)
        dx = f/df
        rtnewtonOLDTMP = rtnewtonOLD - dx
        write(*,*) "rtnewtonOLD, j:", rtnewtonOLD, j
        write(*,*) "rtnewtonOLDTMP:", rtnewtonOLDTMP
        write(*,*) "     df, dx:", df, dx

*       Check that the new value does not suggest a negative value        
 50     IF (rtnewtonOLDTMP .LT. 0.0d00) THEN
          write(*,*) "ADJUSTING for -ive value!"
          
          rtnewtonOLD = rtnewtonOLD / 2.0d00
          
          CALL funcd(rtnewtonOLD, f, df)
          dx = f/df
          rtnewtonOLDTMP = rtnewtonOLD - dx

          GOTO 50
        ENDIF
        
        rtnewtonOLD = rtnewtonOLDTMP
        write(*,*) "rtnewtonOLD, j:", rtnewtonOLD, j
        write(*,*) "     df, dx:", df, dx
        
*        IF ( (x1 - rtnewtonOLD) * (rtnewtonOLD - x2) .LT. 0.0d00) THEN
*           rtnewtonOLD = DABS(rtnewtonOLD)
*        ENDIF
        
        IF ( DABS(dx) .LT. xacc) RETURN 
*       Convergence!
      ENDDO

      RETURN
      END