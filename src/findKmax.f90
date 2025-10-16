SUBROUTINE findKmax(i, kmax, tmax, mmax, mfirst, startPoint) 
  USE tweedie_params_mod

  IMPLICIT NONE
  
  ! Arguments
  REAL(KIND=8), INTENT(OUT)     :: kmax, tmax
  REAL(KIND=8), INTENT(IN)      :: startPoint
  INTEGER, INTENT(OUT)          :: mmax, mfirst
  INTEGER, INTENT(IN)           :: i
  
  
  ! --- CRITICAL FIXES: INTERFACES FOR EXTERNAL BIND(C) ROUTINES ---
  
  INTERFACE
      ! A) Define the REQUIRED 4-argument signature for the function pointer (funcd)
      SUBROUTINE funcd_signature(i_in, t, f, df)

          IMPLICIT NONE
          INTEGER, INTENT(IN)       :: i_in
          REAL(KIND=8), INTENT(IN)  :: t
          REAL(KIND=8), INTENT(OUT) :: f, df
      END SUBROUTINE funcd_signature

      ! B) 1. The root solver (rtsafe)
      SUBROUTINE rtsafe(i_in, funcd, x1, x2, xstart, xacc, root) 

          ! FIX 1 & 2: Declare funcd using the signature, remove POINTER
          PROCEDURE(funcd_signature) :: funcd 
          
          INTEGER, INTENT(IN)       :: i_in
          REAL(KIND=8), INTENT(IN)  :: x1, x2, xstart, xacc
          REAL(KIND=8), INTENT(OUT) :: root
      END SUBROUTINE rtsafe
      
      ! C) Define the rtnewton solver (used on line 99)
      SUBROUTINE rtnewton(i_in, funcd, x1, x2, xstart, xacc, root)

          PROCEDURE(funcd_signature) :: funcd 
          INTEGER, INTENT(IN)       :: i_in
          REAL(KIND=8), INTENT(IN)  :: x1, x2, xstart, xacc
          REAL(KIND=8), INTENT(OUT) :: root
      END SUBROUTINE rtnewton
      
      ! D) Define findImdkZero (the actual 4-argument routine being passed)
      SUBROUTINE findImdkZero(i_in, t, f, df)

          IMPLICIT NONE
          INTEGER, INTENT(IN)       :: i_in
          REAL(KIND=8), INTENT(IN)  :: t
          REAL(KIND=8), INTENT(OUT) :: f, df
      END SUBROUTINE findImdkZero
      
      SUBROUTINE findImk(i_in, t_in, kmax_out) 

          IMPLICIT NONE
          INTEGER, INTENT(IN)         :: i_in
          REAL(KIND=8), INTENT(IN)    :: t_in
          REAL(KIND=8), INTENT(OUT)   :: kmax_out
      END SUBROUTINE findImk
      
      ! F) Define myfloor (used on line 72, 107)
      INTEGER FUNCTION myfloor(x) 

          IMPLICIT NONE
          REAL(KIND=8), INTENT(IN) :: x
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
      
      ! FIX 4: The call site is now correct because funcd's type is explicitly defined.
      CALL rtsafe(i, findImdkZero, kmaxL, kmaxR, startPoint, aimrerr, tmax)
      
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

      CALL rtnewton(i, findImdkZero, 0.00d0, startPoint * 30.0d0, startPoint, aimrerr, tmax)
      ! funcd returns the fn value, and derivative value

      ! Find kmax, mmax
      CALL findImk(i, tmax, kmax)
      mmax = myfloor(kmax/pi)
      mfirst = mmax
      
    END IF
  END IF

  RETURN

END SUBROUTINE findKmax
