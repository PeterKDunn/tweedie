SUBROUTINE findKmax(i, kmax, tmax, mmax, mfirst, startPoint)  BIND(C, NAME='findKmax')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Arguments
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: kmax, tmax
  REAL(KIND=C_DOUBLE), INTENT(IN)     :: startPoint
  INTEGER(C_INT), INTENT(OUT)         :: mmax, mfirst
  INTEGER(C_INT), INTENT(IN)          :: i
  
  ! Local Variables
  REAL(KIND=8)        :: pi
  
 ! Fix: rtnewton and rtsafe are function calls that return a C_DOUBLE
  REAL(KIND=C_DOUBLE)        :: rtnewton, rtsafe
  REAL(KIND=C_DOUBLE)        :: findImdkZero, aimrerr
  
  ! Variables used internally
  REAL(KIND=8)   :: kmaxL, kmaxR
  REAL(KIND=8)   :: current_y, current_mu, current_phi
  INTEGER        :: myfloor

  EXTERNAL findImdkZero, myfloor, rtnewton, rtsafe

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  ! startPoint = 0.0d0
  
  aimrerr = 1.0d-09
  pi = 4.0D0 * DATAN(1.0D0)

  IF ( CpSmall) THEN
    IF (current_y .GE. current_mu) THEN
      ! Cy >= Cmu and 1 < p < 2
      mmax = 0
      tmax = 0.0d00
      kmax = 0.0d00
      ! startPoint = pi / current_y + 0.05d00
    ELSE
      ! Cy < Cmu and 1 < p < 2
  WRITE(*,*) "IN findKmax: startPoint", startPoint
  WRITE(*,*) "About to call rtnewton/rtsafe"
      tmax = rtsafe(findImdkZero, kmaxL, kmaxR, startPoint, aimrerr)
      ! funcd returns the fn value, and derivative value

      CALL findImk(i, tmax, kmax)
      mmax = myfloor(kmax/pi)
      mfirst = mmax
    END IF
  ELSE
    ! IF p > 2
    IF (current_y .GE. current_mu) THEN
      ! Cy >= Cmu and p > 2

      mmax = 0
      mfirst = -1
      tmax = 0.0d00
      kmax = 0.0d00
    ELSE
      ! Cy < Cmu and p > 2

      tmax = rtnewton(findImdkZero,      &
                      0.00d00,           &
                      startPoint * 30.0d00, &
                      startPoint,        &
                      aimrerr)
      ! funcd returns the fn value, and derivative value

      ! Find kmax, mmax
      CALL findImk(i, tmax, kmax)
      mmax = myfloor(kmax/pi)
      mfirst = mmax
      
    END IF
  END IF

  RETURN

END SUBROUTINE findKmax


