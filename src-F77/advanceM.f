      SUBROUTINE advanceM(mmax, m, mOld, leftOfMax, flip)
    
      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      INTEGER mOld, m, mmax
      LOGICAL leftOfMax, pSmall, flip

      mOld = m
      flip = .FALSE.
      
      IF (pSmall) THEN
*       IS THIS a separate IF??? OR use same one???

        IF (Cy .GE. Cmu) THEN
*         Always heading downwards, so easy
          m = m - 1
        ELSE
          IF (leftOfMax) THEN
            IF (m .EQ. mmax) THEN
*             Move to the other side
              leftOfMax = .FALSE.
              flip = .TRUE.
            ELSE
              mOld = m
              m = m + 1
              leftOfMax = .TRUE.
            ENDIF
          ELSE
*           Can always just reduce m by one when to the the RIGHT of the maximum
             m = m - 1
             leftOfMax = .FALSE.
          ENDIF
        
        ENDIF


      ELSE
*       FOR: p is large
        IF (Cy .GE. Cmu) THEN
*         Always heading downwards, so easy
          m = m - 1
        ELSE
          IF (leftOfMax) THEN
            IF (m .EQ. mmax) THEN
*             Move to the other side
              leftOfMax = .FALSE.
              flip = .TRUE.
            ELSE
              mOld = m
              m = m + 1
              leftOfMax = .TRUE.
            ENDIF
          ELSE
*           Can always just reduce m by one when to the the RIGHT of the maximum
             m = m - 1
             leftOfMax = .FALSE.
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      
      