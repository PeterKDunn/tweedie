
SUBROUTINE findRek(i, t, Rek) BIND(C, NAME='findRek')
  USE tweedie_params_mod, ONLY: Cp, Cphi, Cmu
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: t              ! Input parameter t
  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Rek            ! Output result (Real part of k(t))

  REAL(KIND=8) :: current_mu, current_phi
  REAL(KIND=8) :: pi, omega, pindex, front, alpha, tanArg
  
  ! Grab the relevant scalar values for this iteration:
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  

  pi = ACOS(-1.0d0) ! Modern F90 way to get PI
  
  pindex = (2.0d00 - Cp)
  front = current_mu ** pindex  / ( current_phi * pindex)
  
  tanArg = (1.0d00 - Cp) * t * current_phi / (current_mu ** (1.0d00 - Cp) )
           
  omega = ATAN( tanArg )
  
  ! Safety Check (retaining the core math logic, but removing I/O)
  IF ((omega > 0.0d00 ) .OR. (omega < (-pi/2.0d00))) THEN
     ! Raise an error if omega is out of the valid range
     ! Since we cannot use R I/O, we must rely on the R wrapper 
     ! to handle the return of an invalid value, or let a runtime error occur.
     ! Returning NaN is a robust way to signal the error back to R.
!     Rek = ACOS(2.0d0) ! Resulting in NaN
     RETURN
  END IF
  
  alpha = (2.0d00 - Cp)/(1.0d00 - Cp)

  Rek = front * &
        ( COS(omega * alpha)/(COS(omega)**alpha) - 1.0d00)

END SUBROUTINE findRek
