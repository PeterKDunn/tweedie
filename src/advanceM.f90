SUBROUTINE advanceM(i, m, mmax, mOld, leftOfMax, flip)

  USE tweedie_params_mod, ONLY: Cy, Cmu, Cphi
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  INTEGER(C_INT), INTENT(IN)    :: mmax, i    ! Maximum value m can take, and current index
  INTEGER(C_INT), INTENT(INOUT) :: m          ! M index (used for calculation and C-binding)
  INTEGER(C_INT), INTENT(OUT)   :: mOld       ! Previous index value
  INTEGER(C_INT), INTENT(INOUT) :: leftOfMax  ! True if on the left side of kmax
  INTEGER(C_INT), INTENT(INOUT) :: flip       ! True if cross from left to right
  
  REAL(KIND=C_DOUBLE)                  :: current_y, current_mu, current_phi


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  mOld = m
  flip = 0
  
  IF (current_y .GE. current_mu) THEN
    ! Always heading downwards (away from kmax), so easy
    m = m - 1 
  ELSE
    ! We have a maximum (kmax) to consider
    IF (leftOfMax .EQ. 1) THEN
      IF (m == mmax) THEN 
        ! Move to the other side of the maximum
        leftOfMax = 0
        flip = 1
      ELSE
        ! Continue towards the maximum
        m = m + 1 
        ! mOld is already saved before the IF block
        leftOfMax = 1 ! Still on the left side
      END IF
    ELSE
      ! When on the RIGHT of the maximum, can always just reduce m by one
      m = m - 1 
      leftOfMax = 0
    END IF
    
  END IF

END SUBROUTINE advanceM
