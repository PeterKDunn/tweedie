      SUBROUTINE findExactZeros(zeroL, zeroR, zeroSP, zero)
*     Finds the exact zeros, to integrate between

      IMPLICIT NONE
      DOUBLE PRECISION zeroR, zeroL, zero, zeroSP
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      DOUBLE PRECISION xacc, rtnewton
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      EXTERNAL findImkM, findZeroSmallp
      
      xacc = 1.0d-12
*     Sometimes, the curve is quite flat and accuracy is hard to 
*     achieve, but not vital anyway

      IF (pSmall) THEN
        zero = rtnewton(findZeroSmallp, zeroL, zeroR, zeroSP, xacc)
      ELSE
        zero = rtnewton(findImkM, zeroL, zeroR, zeroSP, xacc)
      ENDIF
      
      RETURN
      END
