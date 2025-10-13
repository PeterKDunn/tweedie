
SUBROUTINE findImkdd(i, t, Imkdd) BIND(C, NAME='findImkdd')
  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Arguments
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkdd
  
  ! Local Variables
  REAL(KIND=8)    :: front, omega, pindex
  REAL(KIND=8)    :: current_mu, current_phi

  ! Grab the relevant scalar values for this iteration:
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  ! F90 constants use the _8 suffix for double precision
  pindex = Cp / (1.0_8 - Cp)
  front = -current_phi * current_mu ** (Cp/(1.0_8 - Cp))
  
  ! The calculation using F90 syntax
  ! DATAN is replaced by the standard Fortran D-prefix function, DATAN(X) 
  ! or simply ATAN(X) with REAL(KIND=8) arguments. We use DATAN for clarity.
  omega = DATAN( ( (1.0_8 - Cp) * t * current_phi) / (current_mu ** (1.0_8 - Cp) ) )

  Imkdd = front * (DSIN(omega * pindex) / (DCOS(omega) ** pindex) )

END SUBROUTINE findImkdd
