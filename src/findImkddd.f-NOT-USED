      SUBROUTINE findImkddd(t, Imdddk)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy
      DOUBLE PRECISION Imdddk, front, alpha, omega, pindex, t
      COMMON /params/ Cp, Cy, Cmu, Cphi

      write(*,*) "- Finding Im(k''(t))", Cy

      front = -1.0d00 * Cphi ** 2.0d00 * Cp
      omega = DATAN( ( (1.0d00 - Cp) * t * Cphi)/
     &               (Cmu ** (1.0d00 - Cp) ) )
      alpha = (2.0d00 - Cp)/(1.0d00 - Cp)
      pindex = (2.0d00 * Cp - 1.0d00) / ( 1.0d00 - Cp)

      Imdddk = front *
     &        DCOS( omega * pindex) /
     &        (DCOS(omega) ** pindex)

      RETURN 
      END