SUBROUTINE findKmax(i, kmax, tmax, mmax, mfirst, leftOfMax) 
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Arguments
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: kmax, tmax
  INTEGER(C_INT), INTENT(OUT)         :: mmax, mfirst, leftOfMax
  INTEGER(C_INT), INTENT(IN)          :: i

  
  INTERFACE
    FUNCTION findKmaxSP(j) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE  
      REAL(KIND=C_DOUBLE)          :: findKmaxSP
      INTEGER, INTENT(IN)   :: j
    END FUNCTION findKmaxSP


    SUBROUTINE improveKZeroBounds(i, m, leftOfMax, startx, xL, xR)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), INTENT(IN)       :: i, m, leftOfMax
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: startx
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: xL, xR
    END SUBROUTINE improveKZeroBounds
    
    
    SUBROUTINE funcd_signature(i_in, t, f, df) BIND(C)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN) :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
    END SUBROUTINE funcd_signature


    SUBROUTINE rtsafe(i_in, funcd, xstart, x1, x2, xacc, root) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      PROCEDURE(funcd_signature) :: funcd 
    
      INTEGER(C_INT), INTENT(IN) :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN) :: x1, x2, xstart, xacc
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: root
    END SUBROUTINE rtsafe
      

    SUBROUTINE rtnewton(i_in, funcd, xstart, x1, x2, xacc, root) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      PROCEDURE(funcd_signature) :: funcd 
      INTEGER(C_INT), INTENT(IN) :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN) :: x1, x2, xstart, xacc
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: root
    END SUBROUTINE rtnewton
      

    SUBROUTINE findImkdZero(i_in, t, f, df) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN) :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
    END SUBROUTINE findImkdZero
      

    SUBROUTINE findImk(i_in, t_in, kmax_out) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN) :: t_in
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: kmax_out
    END SUBROUTINE findImk
      

    INTEGER(C_INT) FUNCTION myfloor(x) 
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(IN) :: x
    END FUNCTION myfloor
      
  END INTERFACE


  ! Local Variables
  REAL(KIND=C_DOUBLE)     :: pi, omega_Inflection, t_SP
  REAL(KIND=C_DOUBLE)     :: aimrerr, kmaxL, kmaxR, omega_SP, ratio
  REAL(KIND=C_DOUBLE)     :: current_y, current_mu, current_phi
  
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! --- Initialization ---
  aimrerr = 1.0E-09_C_DOUBLE
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

  ! Find starting points
  IF (current_y .GE. current_mu) THEN
    ! CASE: IF y >= mu, and Im k(t) heads down
    ! Nothing to do; easy-peasy:
    mmax = 0
    mfirst = -1
    kmax = 0.0_C_DOUBLE
    tmax = 0.0_C_DOUBLE
    leftOfMax = 0
  ELSE
    ! CASE: IF y < mu: trickier, esp. with 1 < p < 2
    ! Good starting point often needed
!    IF (Cpsmall) THEN
!      t_SP = findKmaxSP(i)
!      CALL improveKmaxSPbounds(i, t_SP, kmaxL, kmaxR)
!    ELSE
      ratio = current_y / current_mu
      omega_SP = -1
      IF (ratio .LT. 0.1_C_DOUBLE) omega_SP = -1.5_C_DOUBLE
      IF (ratio .GT. 0.9_C_DOUBLE) omega_SP = 0.01_C_DOUBLE
      
      t_SP = current_mu ** (1.0_C_DOUBLE - Cp) * DTAN(omega_SP) / ( ( 1.0_C_DOUBLE - Cp) * current_phi)
      kmaxL = t_SP * 0.5_C_DOUBLE
      kmaxR = t_SP * 2.0_C_DOUBLE
!    END IF
    ! Now find kmax and friends
    CALL rtnewton(i,              &
                  findImkdZero,   &
                  t_SP,           &
                  kmaxL,          &
                  kmaxR,          &
                  aimrerr,        &
                  tmax)

    ! Find kmax, mmax
    CALL findImk(i, tmax, kmax)
    mmax = myfloor(kmax/pi)
    mfirst = mmax
    
    leftOfMax = .TRUE.
    IF (mmax .EQ. 0) THEN
      leftOfMax = .FALSE.
    END IF
    
    END IF

  RETURN

END SUBROUTINE findKmax
