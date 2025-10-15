
SUBROUTINE findImk(i, t, Imk) 
      
  USE tweedie_params_mod

  IMPLICIT NONE
      
  ! Arguments
  REAL(KIND=8), INTENT(IN)     :: t
  REAL(KIND=8), INTENT(OUT)    :: Imk
  INTEGER, INTENT(IN)          :: i
  
  REAL(KIND=8) :: tanArg, omega, front, alpha, pi
  REAL(KIND=8) :: current_y, current_mu, current_phi
      
  pi = 4.0D0 * DATAN(1.0D0)
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  front = current_mu ** (2.0d00 - Cp) / ( current_phi * (2.0d00 - Cp))
  tanArg = (1.0d00 - Cp) * t * current_phi  / (current_mu ** (1.0d00 - Cp) )
  omega = DATAN( tanArg )

  IF ((omega .GT. 0.0d00 ) .OR.    &    
      (omega .LT. (-pi/2.0d00)) ) THEN
         write(*,*) "ERROR (FindImk): Omega out of range:", omega
         write(*,*) "Argument: ", tanArg
         write(*,*) "       t: ", t
         STOP
  END IF
  alpha = (2.0d00 - Cp)/(1.0d00 - Cp)

  Imk = front *   &
        DSIN(omega * alpha)/DCOS(omega) ** alpha - t * current_y

  RETURN

END SUBROUTINE findImk
