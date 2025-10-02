      SUBROUTINE findLambda(lambda)

      IMPLICIT NONE
      DOUBLE PRECISION lambda, Cp, Cmu, Cphi, Cy
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      lambda = 0.0d0
      IF (pSmall) THEN
        lambda = -1.00d00 * (Cmu ** (2.0d00 - Cp) ) /
     &                      (Cphi * (1.0d00 - Cp) )
      ENDIF
      
      RETURN
      END