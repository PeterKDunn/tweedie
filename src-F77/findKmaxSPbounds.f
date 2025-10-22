
      SUBROUTINE improveKmaxSPBounds(startTKmax, kmaxL, kmaxR)

*     For the case 1 < p < 2, and y < mu, finding a good
*     starting points and good bounds is often crucial
*     for finding kmax.
*     So here we try to find better bounds, by checking 
*     the slope of Im k(t).

      IMPLICIT NONE
      DOUBLE PRECISION boundL, boundR, slope
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      DOUBLE PRECISION kmaxL, kmaxR, SPslope, StartTKmax
      DOUBLE PRECISION oldBoundL, oldBoundR
      INTEGER maxSearch
      LOGICAL pSmall, keepSearching
      EXTERNAL findImkdZero, myfloor, rtnewton, rtsafe
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      
*     FIND the slope of the starting point
*     Positive: to the *left* of the max.
*     Negative: to the *right* of the max.

      CALL findImkd(startTKmax, SPslope)

*     ************** LOWER BOUND
*     If slope at SP is *positive*, only need to creep to the right      

      boundL = startTKmax
      IF (SPslope .LE. 0.0d0) THEN
        keepSearching = .TRUE.
  88    IF (keepSearching) THEN
*        If slope at SP is negative, take bold steps left to find lower bound 
          boundL = boundL / 2.0d0
*         NOTE: this will never go negative

          CALL findImkd(boundL, slope)

          IF (slope .GT. 0.0d0 ) THEN
*           Found a lower bound
            keepSearching = .FALSE.
          ENDIF 
          GOTO 88
        ENDIF
      ENDIF

*     Now creep to the right from  boundL   
      keepSearching = .TRUE.
      
 55   IF (keepSearching) THEN
*       If slope at SP is positive, creep to right to improve lower bound 
        oldBoundL = boundL
        boundL = boundL * 1.10d0
        CALL findImkd(boundL, slope)

        IF (slope .LT. 0.0d0 ) THEN
*         Gone too far! Keep previous bound
          keepSearching = .FALSE.
          boundL = oldBoundL
        ENDIF
        GOTO 55
      ENDIF
      kmaxL = boundL
*      write(*,*) "LBOUND: ", kmaxL

*     ************** UPPER BOUND
      CALL findImkd(startTKmax, SPslope)
      boundR = startTKmax
*      write(*,*) "  - rbound (start): ", boundR
      
*     If slope at SP is *negative*, only need to creep to the left      
      IF (SPslope .GT. 0.0d0) THEN
        boundR = startTKmax
        keepSearching = .TRUE.
        
  28    IF (keepSearching) THEN
*        If slope at SP is positive, take bold steps right to find lower bound 
          boundR = boundR * 1.5d0
*      write(*,*) "  - rbound (bold right): ", boundR

          CALL findImkd(boundR, slope)

          IF (slope .LT. 0.0d0 ) THEN
*           Found an upper bound
            keepSearching = .FALSE.
          ENDIF 
          GOTO 28
        ENDIF
      ENDIF

*     Now creep to the left from  boundR   
      keepSearching = .TRUE.
      
 65   IF (keepSearching) THEN
*       If slope at SP is negative, creep to left to improve upper bound 
        oldBoundR = boundR 
        boundR = boundR * 0.90d0
*      write(*,*) "  - rbound (creep left): ", boundR
*       NOTE: this will never go negative
        CALL findImkd(boundR, slope)

        IF (slope .GT. 0.0d0 ) THEN
*         Gone too far! Keep previous bound
          keepSearching = .FALSE.
          boundR = oldBoundR
        ENDIF
        GOTO 65
      ENDIF
      kmaxR = boundR 
*      write(*,*) "RBOUND: ", kmaxR

      RETURN
      END
