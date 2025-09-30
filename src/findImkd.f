      SUBROUTINE findImkd(t, Imdk)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, t
      DOUBLE PRECISION Imdk, omega, pindex
      COMMON /params/ Cp, Cy, Cmu, Cphi

      pindex = 1.0d00 / (1.0d00 - Cp)
      omega = DATAN( ( (1.0d00 - Cp) * t * Cphi)/
     &              (Cmu ** (1.0d00 - Cp) ) )
     
      Imdk = Cmu *
     &       DCOS( omega * pindex ) /
     &       (DCOS(omega)**pindex) - Cy

      RETURN
      END
      
