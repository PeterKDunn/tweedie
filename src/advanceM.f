      SUBROUTINE advanceM(mmax, m, mOld, leftOfMax)
    
      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy
      COMMON /params/ Cp, Cy, Cmu, Cphi
      INTEGER mOld, m, mmax
      LOGICAL leftOfMax

      IF (Cp. GT. 2.0d00 ) THEN
        IF (Cy .GE. Cmu) THEN
*         Always heading downwards, so easy
          m = m - 1
        ELSE
          IF (leftOfMax) THEN

            m = m + 1
            IF (m .EQ. mmax ) THEN
              leftOfMax = .FALSE.
            ENDIF
          ELSE
            m = m - 1
          ENDIF
        ENDIF
      ELSE
*       IS THIS a separate IF???
      ENDIF
      
      RETURN
      END
      
      