      SUBROUTINE findKmax(kmax, tmax, mmax, mfirst, startPoint)

      IMPLICIT NONE
      DOUBLE PRECISION kmax, tmax, startPoint, pi, rtnewton
      DOUBLE PRECISION findImdkZero, aimrerr
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi, myfloor
      INTEGER mmax, mfirst
      LOGICAL pSmall
      EXTERNAL findImdkZero, myfloor
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      
      aimrerr = 1.0d-09
      pi = 4.0d00 * DATAN(1.0d00)

      IF ( pSmall) THEN
        IF (Cy .GE. Cmu) THEN
*         Cy >= Cmu and 1 < p < 2
          mmax = 0
          tmax = 0.0d00
          kmax = 0.0d00
          startPoint = pi / Cy + 0.25d00
        ELSE
*         Cy < Cmu and 1 < p < 2
          tmax = rtnewton(findImdkZero, 
     &                  0.0d00, 
     &                  startPoint * 2.0d00, 
     &                  startPoint, aimrerr)
*         funcd returns the fn value, and derivative value
          mmax = 0
          tmax = 0.0d00
          kmax = 0.0d00
        ENDIF
      ELSE
              IF (Cy .GE. Cmu) THEN
*         Cy >= Cmu and p > 2

          mmax = 0
          mfirst = -1
          tmax = 0.0d00
          kmax = 0.0d00
        ELSE
*         Cy < Cmu and p > 2

*      write(*,*) "*** (fKM) l=0, right=", startPoint * 30.0d0
          tmax = rtnewton(findImdkZero, 
     &                    0.00d00, 
     &                    startPoint * 30.0d00, 
     &                    startPoint, 
     &                    aimrerr)
*         funcd returns the fn value, and derivative value

*         Find kmax, mmax
          CALL findImk(tmax, kmax)
          mmax = myfloor(kmax/pi)
          mfirst = mmax
          
        ENDIF
      ENDIF

      RETURN
      END

