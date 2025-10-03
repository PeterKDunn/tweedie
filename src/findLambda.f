      SUBROUTINE findLambda(lambda)

      IMPLICIT NONE
      DOUBLE PRECISION lambda, Cp, Cmu, Cphi, Cy
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      lambda = 0.0d0
      IF (pSmall) THEN
        lambda = (Cmu ** (2.0d00 - Cp) ) /
     &             (Cphi * (2.0d00 - Cp) )
*       NOTE: No negative sign in front: P(Y=0) = exp(-lambda)
      ENDIF
      
      RETURN
      END