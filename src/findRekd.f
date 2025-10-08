      SUBROUTINE findRekd(t, Redk)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, t
      DOUBLE PRECISION Redk, omega, pindex
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      pindex = 1.0d00 / (1.0d00 - Cp)
      omega = DATAN( ( (1.0d00 - Cp) * t * Cphi)/
     &               (Cmu ** (1.0d00 - Cp) ) )
     
*      write(*,*) "  - findImkd on entry:"
*      write(*,*) "  -", Cp, Cy, Cmu, Cphi
*      write(*,*) "  AND at t = ", t
*      write(*,*) "  -", pindex, omega
      Redk = Cmu *
     &       DSIN( omega * pindex ) /
     &       (DCOS(omega)**pindex) 
*      write(*,*) "  -  findRekd on exit:", Imdk

      RETURN
      END
      
