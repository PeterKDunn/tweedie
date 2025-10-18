SUBROUTINE findRek(i, t, Rek)
  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(IN)    :: t
  INTEGER(C_INT), INTENT(IN)         :: i
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: Rek

  REAL(KIND=C_DOUBLE) :: current_mu, current_phi
  REAL(KIND=C_DOUBLE) :: pi, omega, pindex, front, alpha, tanArg
  
  ! Grab the relevant scalar values for this iteration:
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  

  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  pindex = (2.0d00 - Cp)
  front = current_mu ** pindex  / ( current_phi * pindex)
  tanArg = (1.0d00 - Cp) * t * current_phi / (current_mu ** (1.0d00 - Cp) )
  omega = DATAN( tanArg )
  
  ! Safety check
  IF ((omega > 0.0d00 ) .OR. (omega < (-pi/2.0d00))) THEN
     ! Error!
     WRITE(*,*) "Error: omga out of bounds"
     RETURN
  END IF
  
  alpha = (2.0d00 - Cp)/(1.0d00 - Cp)
  Rek = front * ( DCOS(omega * alpha)/(DCOS(omega)**alpha) - 1.0d00 )

END SUBROUTINE findRek
