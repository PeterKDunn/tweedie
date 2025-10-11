      SUBROUTINE findImk(t, Imk)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, pi, tanArg
      DOUBLE PRECISION t, Imk, omega, front, alpha
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      pi = 4.0d00 * DATAN(1.0d0)

      front = Cmu ** (2.0d00 - Cp) / ( Cphi * (2.0d00 - Cp))
      tanArg = (1.0d00 - Cp) * t * Cphi/
     &               (Cmu ** (1.0d00 - Cp) )
      omega = DATAN( tanArg )

      IF ((omega .GT. 0.0d00 ). OR. 
     &    (omega .LT. (-pi/2.0d00)) ) THEN
         write(*,*) "ERROR (FindImk): Omega out of range:", omega
         write(*,*) "Argument: ", tanArg
         write(*,*) "    t: ", t
         STOP
      ENDIF
      alpha = (2.0d00 - Cp)/(1.0d00 - Cp)

      Imk = front *
     &      DSIN(omega * alpha)/DCOS(omega) ** alpha -
     &      t * Cy

      RETURN 
      END