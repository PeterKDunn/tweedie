
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


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  front = current_mu ** (2.0_C_DOUBLE - Cp) / ( current_phi * (2.0_C_DOUBLE - Cp))
  tanArg = (1.0_C_DOUBLE - Cp) * t * current_phi  / (current_mu ** (1.0_C_DOUBLE - Cp) )
  omega = DATAN( tanArg )
!write(*,*) "t  (1500)  = ", t
!write(*,*) "FRONT (-0.2357023)  = ", front
!write(*,*) "tanArg (-21213.2) = ", tanArg
!write(*,*) "OMEGA (-1.570749) = ", omega

  IF ((omega .GT. 0.0_C_DOUBLE ) .OR.    &    
      (omega .LT. (-pi/2.0_C_DOUBLE)) ) THEN
    write(*,*) "ERROR (FindImk): Omega out of range:", omega
!    write(*,*) "Argument: ", tanArg
!    write(*,*) "       t: ", t
    STOP
  END IF
  alpha = (2.0_C_DOUBLE - Cp)/(1.0_C_DOUBLE - Cp)

  Imk = front *   &
        DSIN(omega * alpha)/(DCOS(omega) ** alpha) - t * current_y

!write(*,*) "Imk (30.2101) = ", Imk

!IF ( (t .GT. 1450) .AND. (t .LT. 1550 ) ) STOP

RETURN

END SUBROUTINE findImk
