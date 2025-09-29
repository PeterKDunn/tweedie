      FUNCTION rtsafe(funcd, x1, x2, xacc)
*     From Numerical Recipes
*     Find a root between x1 and x2
*     From Numerical Recipes.
*     Newton-Raphson, but with bisection

      INTEGER MAXIT, numberTries
      DOUBLE PRECISION rtsafe, x1, x2, xacc
      EXTERNAL funcd
      PARAMETER (MAXIT=10) 
*     Maximum allowed number of iterations.

      INTEGER j
      DOUBLE PRECISION df, dx, dxold, f, fh, fl, temp, xh, xl
      LOGICAL tryHarder, boundsOK

*     Evaluate function at the endpoints
      CALL funcd(x1, fl, df)
      CALL funcd(x2, fh, df)
      
      boundsOK = .TRUE.
      write(*,*) "  rtsafe: STARTING BOUNDS:", x1, x2
      write(*,*) "  rtsafe: GIVING :", fl, fh
      IF ( (fl .GT. 0.0d00 .AND. fh .GT. 0.0d00) .OR.
     &     (fl .LT. 0.0d00 .AND. fh .LT. 0.0d00) ) THEN
         tryHarder = .TRUE.
         boundsOK = .FALSE.
         write(*,*) "rtsafe: Initial bounds dont bracket zero..."
      ENDIF
      
      numberTries = 0
 19   IF ( tryHarder ) THEN
         numberTries = numberTries + 1
         write(*,*) "  rtsafe: Number of tries:", numberTries
         tryHarder = .FALSE.
         boundsOK = .FALSE.
         write(*,*) "  rtsafe: ERROR: root not bracketed:", fl, fh
         x1 = MAX(0.0d00, x1 * 0.75d00)
         x2 = x2 * 1.25d00
        CALL funcd(x1, fl, df)
        CALL funcd(x2, fh, df)
        
        IF ( (fl .GT. 0.0d00 .AND. fh .GT. 0.0d00) .OR.
     &       (fl .LT. 0.0d00 .AND. fh .LT. 0.0d00) ) THEN
        write(*,*) "   rtsafe: NOPE!"
           tryHarder = .TRUE.
           boundsOK = .FALSE.
        ELSE
           write(*,*) "   rtsafe: >> Success bounding the zero!"
           tryHarder = .FALSE.
           boundsOK = .TRUE.
        ENDIF
        IF ( (ABS(x1) .GT. 1.0d25) .OR. (ABS(x2) .GT. 1.0d25 ) ) THEN
           boundsOK = .FALSE.
           tryHarder = .FALSE.
        ENDIF
        
        IF (numberTries .GT. MAXIT) THEN
           boundsOK = .FALSE.
           tryHarder = .FALSE.
        ENDIF
        
        write(*,*) "  rtsafe: RETRY: update: x1=", x1," x2=", x2
        write(*,*) "                 with fl=", fl," fh=", fh
        goto 19 
      ENDIF
      
      IF (.NOT.(boundsOK)) THEN
        write(*,*) "  rtsafe: ERROR: NO GOOD BOUNDS FIND: ERROR"
        write(*,*) "rtsafe:   Stopped with", x1, x2
        write(*,*) "          having", fl, fh
        STOP
      ENDIF

      write(*,*) "  rtsafe: FOUND BOUNDS: x1=", x1, " x2=",x2

      IF (fl .EQ. 0.0d00) THEN
        rtsafe = x1
        RETURN
      ELSE IF (fh .EQ. 0.0d00) THEN
        rtsafe = x2
        RETURN
      ELSE IF (fl .LT. 0.0d00) THEN
*       Orient the search so that f(xl) < 0.
        xl = x1
        xh = x2
      ELSE
        xh = x1
        xl = x2
      ENDIF

*     Initial guess for root
      rtsafe = 0.5d00 * (x1 + x2)

      dxold = abs(x2 - x1)
      dx    = dxold

      CALL funcd(rtsafe, f, df)

      DO j = 1, MAXIT
*       Loop over allowed iterations.
         write(*,*) "------ rtsafe: ITERATING ------ j=", j
         write(*,*) "  rtsafe: Using : fl=", fl," fh=",fh
         write(*,*) "rtsafe: Using : df", df

        if ( ( ( (rtsafe - xh) * df - f) *
     &         ( (rtsafe - xl) * df - f) .GT. 0.0d00 )
     &      .OR.
     &       abs(2.0d00 * f) .GT. abs(dxold * df) ) then
*         Bisect if Newton out of range, or not decreasing fast enough.
          dxold  = dx
          dx     = 0.5d00 * (xh - xl)
          rtsafe = xl + dx

          IF (xl .EQ. rtsafe) THEN
*           Change in root is negligible.
            return
          ELSE
*           Newton step acceptable; take it.
            dxold  = dx
            dx     = f / df
            temp   = rtsafe
            rtsafe = rtsafe - dx
            IF (temp .EQ. rtsafe) RETURN
          ENDIF

        ELSE
*         Newton step was fine
          dxold  = dx
          dx     = f / df
          temp   = rtsafe
          rtsafe = rtsafe - dx
          IF (temp .EQ. rtsafe) RETURN
        ENDIF

*       Convergence test
        IF (abs(dx) .LT. xacc) RETURN

*       Evaluate function at new point
        CALL funcd(rtsafe, f, df)

*       Maintain the bracket
        IF (f .LT. 0.0d00) THEN
          xl = rtsafe
        ELSE
          xh = rtsafe
        ENDIF

      ENDDO

*     If we get here, we exceeded maximum iterations
      write(*,*) ' mrtsafe: exceeding maximum iterations'

      RETURN
      END


*****************************************************************
*****************************************************************