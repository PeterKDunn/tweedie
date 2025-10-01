      SUBROUTINE findImkdd(t, Imddk)
      
      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy
      DOUBLE PRECISION Imddk, front, omega, pindex, t
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      pindex = Cp / (1.0d00 - Cp)
      front = -1.0d00 * Cphi * Cmu**Cp
***** CHANGED FRONT FROM THESES!

      omega = DATAN( ( (1.0d00 - Cp) * t * Cphi)/
     &               (Cmu ** (1.0d00 - Cp) ) )

      Imddk = front *
     &        DSIN( omega * pindex) /
     &        (DCOS(omega) ** (pindex) )

      RETURN
      END