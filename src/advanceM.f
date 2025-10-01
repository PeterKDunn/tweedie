      SUBROUTINE advanceM(mmax, m, mOld, leftOfMax)
    
      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      INTEGER mOld, m, mmax
      LOGICAL leftOfMax, pSmall

      mOld = m
      
      IF (pSmall) THEN

*       IS THIS a separate IF??? OR use same one???

      ELSE
        IF (Cy .GE. Cmu) THEN
*         Always heading downwards, so easy
          m = m - 1
        ELSE
          IF (leftOfMax) THEN
            m = m + 1
            IF ( (m .EQ. mmax ) .AND.
     &           (mOld .EQ. mmax ) ) THEN
              leftOfMax = .FALSE.
            ENDIF
          ELSE
            m = m - 1
          ENDIF
        ENDIF
      ENDIF
      
      RETURN
      END
      
      