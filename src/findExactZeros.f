      SUBROUTINE findExactZeros(zeroL, zeroR, zeroSP, zero)
*     Finds the exact zeros, to integrate between

      IMPLICIT NONE
      INTEGER m
      DOUBLE PRECISION zeroR, zeroL, zero, zeroSP
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      DOUBLE PRECISION xacc, rtnewton
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      COMMON /mparam/ m
      EXTERNAL findImkM
      
      xacc = 1.0d-09
*     Sometimes, the curve is quite flat and accuracy is hard to 
*     achieve, but not vital anyway
      zero = rtnewton(findImkM, zeroL, zeroR, zeroSP, xacc)

      RETURN
      END
