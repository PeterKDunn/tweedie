      SUBROUTINE findImkdd(t, Imkdd)
      
      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy
      DOUBLE PRECISION Imkdd, front, omega, pindex, t
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      pindex = Cp / (1.0d00 - Cp)
      front = -Cphi * Cmu ** (Cp/(1.0d0 - Cp))
***** CHANGED FRONT FROM THESES!

      omega = DATAN( ( (1.0d00 - Cp) * t * Cphi)/
     &               (Cmu ** (1.0d00 - Cp) ) )

      Imkdd = front *
     &        DSIN( omega * pindex) /
     &        (DCOS(omega) ** (pindex) )

      RETURN
      END
      