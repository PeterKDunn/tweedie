!!!! NEED T BRANCH BY p???



SUBROUTINE advanceM(i, m, mmax, mOld, leftOfMax, flip)
  ! --- FIX: Import the global module variable 'm' and rename it to 'global_m'
  ! This allows us to access and update the shared state without conflicting
  ! with the local C-bound argument 'm_index'.
  USE tweedie_params_mod, ONLY: CpSmall, Cy, Cmu, Cphi
  USE ISO_C_BINDING, ONLY: C_INT

  IMPLICIT NONE
  
  ! Arguments (Inputs/Outputs)
  INTEGER, INTENT(IN)    :: mmax, i          ! Maximum value m can take, and current index
  INTEGER, INTENT(INOUT) :: m                ! M index (used for calculation and C-binding)
  INTEGER, INTENT(OUT)   :: mOld             ! Previous index value
  INTEGER, INTENT(INOUT) :: leftOfMax        ! True if on the left side of kmax
  INTEGER, INTENT(INOUT) :: flip             ! True if cross from left to right
  
  REAL(KIND=8)                  :: current_y, current_mu, current_phi

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  mOld = m
  flip = 0
  
  IF (CpSmall) THEN
    ! FOR: 1 < p < 2
    
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
    
  ELSE
    ! FOR: p is large (p > 2 or p < 1)
    
    IF (current_y .GE. current_mu) THEN
      ! Always heading downwards, so easy
      m = m - 1 
    ELSE
      ! We have a maximum (kmax) to consider
      IF (leftOfMax .EQ. 1) THEN
        IF (m == mmax) THEN 
          ! Move to the other side
          leftOfMax = 0
          flip = 1
        ELSE
          ! Continue towards the maximum
          m = m + 1 
          ! mOld is already saved before the IF block
          leftOfMax = 1
        END IF
      ELSE
        ! Can always just reduce m by one when to the the RIGHT of the maximum
        m = m - 1 
        leftOfMax = 0
      END IF
    END IF
  END IF

  ! --- CRITICAL STEP: Update the global module variable 'm'
  ! This ensures the shared state is updated after the calculation.
  
END SUBROUTINE advanceM
