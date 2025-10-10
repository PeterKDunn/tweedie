      SUBROUTINE findKmax(kmax, tmax, mmax, mfirst, startPoint,
     &                    kmaxL, kmaxR)

      IMPLICIT NONE
      DOUBLE PRECISION kmax, tmax, startPoint, pi, rtnewton
      DOUBLE PRECISION findImdkZero, aimrerr, rtsafe
      DOUBLE PRECISION kmaxL, kmaxR
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      INTEGER mmax, mfirst, myfloor
      LOGICAL pSmall
      EXTERNAL findImdkZero, myfloor, rtnewton, rtsafe
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
      write(*,*) "IN findKmax: startPoint", startPoint
      write(*,*) "About to call rtnewton/rtsafe"
*          tmax = rtsafe(findImdkZero, 
*     &                    0.0d00, 
*     &                    startPoint * 20.0d00, 
*     &                    startPoint, aimrerr)
          tmax = rtsafe(findImdkZero, 
     &                 kmaxL, kmaxR, startPoint, aimrerr)
*         funcd returns the fn value, and derivative value
*      write(*,*) "BACK IN findKmax"
*         funcd returns the fn value, and derivative value
*         Find kmax, mmax
          CALL findImk(tmax, kmax)
          mmax = myfloor(kmax/pi)
          mfirst = mmax
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

