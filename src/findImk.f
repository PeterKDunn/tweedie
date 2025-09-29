      SUBROUTINE findImk(t, Imk)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, t
      DOUBLE PRECISION Imk, omega, pindex, front, alpha
      COMMON /params/ Cp, Cy, Cmu, Cphi

      front = Cmu ** (2.0d00 - Cp) / ( Cphi * (2.0d00 - Cp))
      omega = DATAN( ( (1.0d00 - Cp) * t * Cphi)/
     &              (Cmu ** (1.0d00 - Cp) ) )
      alpha = (2.0d00 - Cp)/(1.0d00 - Cp)

      Imk = front *
     &      DSIN(omega * alpha)/DCOS(omega) ** alpha -
     &      t * Cy

      RETURN 
      END