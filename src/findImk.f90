
SUBROUTINE findImk(i, t, Imk) 
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
      
  ! Arguments
  REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Imk
  INTEGER(C_INT), INTENT(IN)          :: i
  
  REAL(KIND=C_DOUBLE) :: tanArg, omega, front, alpha, pi
  REAL(KIND=C_DOUBLE) :: current_y, current_mu, current_phi
      
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
        DSIN(omega * alpha)/(DCOS(omega) ** alpha) - t * current_y

  RETURN

END SUBROUTINE findImk
