      SUBROUTINE advanceM(mmax, kmax, leftSide, mOld, mNew)
    
      IMPLICIT NONE
      DOUBLE PRECISION kmax
      INTEGER mOld, mNew, mmax
      LOGICAL leftSide
      
      if (mOld .EQ. mmax) then
        if (leftSide) then
          mNew = mmax
          leftSide = .FALSE.
        else
          mNew = mOld - 1
        endif
      else
        if (leftSide) then
          mNew = mOld + 1
        else
          mNew = mOld - 1
        endif
      endif
      
      RETURN
      END
      
      