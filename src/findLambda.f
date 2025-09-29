      SUBROUTINE findLambda(lambda, p, mu, phi)

      IMPLICIT NONE
      DOUBLE PRECISION lambda, p, mu, phi

      lambda = -1.00d00 * (mu ** (2.0d00 - p) ) /
     &                    ( phi * (1.0d00 - p) )

      RETURN
      END