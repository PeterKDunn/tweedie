      SUBROUTINE findRek(t, Rek)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, t
      DOUBLE PRECISION omega, pindex, front, alpha
      DOUBLE PRECISION pi, Rek, tanArg
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      pi = 4.0d00 * DATAN(1.0d0)

      pindex = (2.0d00 - Cp)
      front = Cmu ** pindex  / ( Cphi * pindex)
      tanArg = (1.0d00 - Cp) * t * Cphi/
     &               (Cmu ** (1.0d00 - Cp) )
      omega = DATAN( tanArg )
      IF ((omega .GT. 0.0d00 ). OR. 
     &    (omega .LT. (-pi/2.0d00)) ) THEN
         write(*,*) "ERROR (FindImk): Omega out of range:", omega
         write(*,*) "Argument: ", tanArg
         write(*,*) "Argument: ", tanArg
         write(*,*) "    p: ", Cp
         write(*,*) "    phi: ", Cphi
         write(*,*) "    mu: ", Cmu
         write(*,*) "    t: ", t
      ENDIF
      alpha = (2.0d00 - Cp)/(1.0d00 - Cp)

      Rek = front *
     &      ( DCOS(omega * alpha)/(DCOS(omega)**alpha) - 1.0d00)

      RETURN
      END