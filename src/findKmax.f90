SUBROUTINE findKmax(i, kmax, tmax, mmax, mfirst, startPoint) 
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Arguments
  REAL(KIND=C_DOUBLE), INTENT(OUT)     :: kmax, tmax
  REAL(KIND=C_DOUBLE), INTENT(IN)      :: startPoint
  INTEGER(C_INT), INTENT(OUT)          :: mmax, mfirst
  INTEGER(C_INT), INTENT(IN)           :: i
  
  
  ! --- CRITICAL FIXES: INTERFACES FOR EXTERNAL BIND(C) ROUTINES ---
  
  INTERFACE
      ! A) Define the REQUIRED 4-argument signature for the function pointer (funcd)
      SUBROUTINE funcd_signature(i_in, t, f, df) BIND(C)
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
          IMPLICIT NONE
          INTEGER(C_INT), INTENT(IN) :: i_in
          REAL(KIND=C_DOUBLE), INTENT(IN) :: t
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
      END SUBROUTINE funcd_signature

      ! B) 1. The root solver (rtsafe)
      SUBROUTINE rtsafe(i_in, funcd, x1, x2, xstart, xacc, root) 
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
          
          ! FIX 1 & 2: Declare funcd using the signature, remove POINTER
          PROCEDURE(funcd_signature) :: funcd 
          
          INTEGER(C_INT), INTENT(IN) :: i_in
          REAL(KIND=C_DOUBLE), INTENT(IN) :: x1, x2, xstart, xacc
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: root
      END SUBROUTINE rtsafe
      
      ! C) Define the rtnewton solver (used on line 99)
      SUBROUTINE rtnewton(i_in, funcd, x1, x2, xstart, xacc, root) 
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
          PROCEDURE(funcd_signature) :: funcd 
          INTEGER(C_INT), INTENT(IN) :: i_in
          REAL(KIND=C_DOUBLE), INTENT(IN) :: x1, x2, xstart, xacc
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: root
      END SUBROUTINE rtnewton
      
      ! D) Define findImkdZero (the actual 4-argument routine being passed)
      SUBROUTINE findImkdZero(i_in, t, f, df) 
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
          IMPLICIT NONE
          INTEGER(C_INT), INTENT(IN) :: i_in
          REAL(KIND=C_DOUBLE), INTENT(IN) :: t
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
      END SUBROUTINE findImkdZero
      
      ! E) Define findImk (used on line 71, 106)
      SUBROUTINE findImk(i_in, t_in, kmax_out) 
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
          IMPLICIT NONE
          INTEGER(C_INT), INTENT(IN) :: i_in
          REAL(KIND=C_DOUBLE), INTENT(IN) :: t_in
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: kmax_out
      END SUBROUTINE findImk
      
      ! F) Define myfloor (used on line 72, 107)
      INTEGER(C_INT) FUNCTION myfloor(x) 
          USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
          IMPLICIT NONE
          REAL(KIND=C_DOUBLE), INTENT(IN) :: x
      END FUNCTION myfloor
      
  END INTERFACE

  ! Local Variables
  REAL(KIND=8)        :: pi
  
  ! Variables used internally
  REAL(KIND=8)        :: kmaxL, kmaxR, aimrerr
  REAL(KIND=8)        :: current_y, current_mu, current_phi
  
  ! Explicitly define local integer variable
  ! INTEGER(C_INT)      :: mmax_local, mfirst_local ! These are OUT arguments now
  
  ! --- Initialization ---
  
  aimrerr = 1.0d-09
  pi = 4.0D0 * DATAN(1.0D0)

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! Default boundaries for rtsafe (assuming the logic expects these values)
  kmaxL = 0.0d00
  kmaxR = pi / 2.0d00  ! A typical starting point for tan(t) root finding
  
  IF ( CpSmall) THEN
    IF (current_y .GE. current_mu) THEN
      ! Cy >= Cmu and 1 < p < 2
      mmax = 0
      mfirst = 0
      tmax = 0.0d00
      kmax = 0.0d00
      
    ELSE
      ! Cy < Cmu and 1 < p < 2
      CALL rtsafe(i,            &
                  findImkdZero, &
                  kmaxL,        &
                  kmaxR,        &
                  startPoint,   &
                  aimrerr,      &
                  tmax)
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
      CALL rtnewton(i,                    &
                    findImkdZero,         &  ! function to find the zeros for
                    0.00d00,              &  ! Lower bound
                    startPoint * 30.0d00, &  ! Upper bound
                    startPoint,           &  ! Starting value
                    aimrerr,              &  ! The relative error aimed for
                    tmax)                    ! The root
    END IF
  END IF
  
  ! Find kmax, mmax
  CALL findImk(i, tmax, kmax)
  mmax = myfloor(kmax/pi)
  mfirst = mmax

  RETURN

END SUBROUTINE findKmax
