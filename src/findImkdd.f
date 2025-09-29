      SUBROUTINE findImkdd(t, Imddk)
      
      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy
      DOUBLE PRECISION Imddk, front, alpha, omega, pindex, t
      COMMON /params/ Cp, Cy, Cmu, Cphi

      pindex = Cp / (1.0d00 - Cp)
      front = -1.0d00 * Cphi * Cmu**pindex
      omega = DATAN( ( (1.0d00 - Cp) * t * Cphi)/
     &               (Cmu ** (1.0d00 - Cp) ) )
      alpha = (2.0d00 - Cp)/(1.0d00 - Cp)

      Imddk = front *
     &        DSIN( omega * pindex) /
     &        (DCOS(omega) ** (pindex) )

      RETURN
      END