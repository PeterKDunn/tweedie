SUBROUTINE improveKZeroBounds(i, m, leftOfMax, startZero, zeroL, zeroR)
  USE tweedie_params_mod, ONLY: Cphi, Cmu, Cy
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Inputs
  REAL(KIND=C_DOUBLE), INTENT(IN)    :: startZero
  INTEGER(C_INT), INTENT(IN)         :: i, m, leftOfMax
  ! Outputs
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: zeroL, zeroR

  ! --- Local Variables ---
  REAL(KIND=C_DOUBLE)     :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE)     :: boundL, boundR, valueL, valueR, SPvalue, multiplier
  REAL(KIND=C_DOUBLE)     :: oldBoundL, oldBoundR, df
  LOGICAL                 :: keepSearching
  
  EXTERNAL findImkM
  

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! Set multipier
  write(*,*) "LEFT OF MAX =", leftOfMax
  IF (leftOfMax .EQ. 1) THEN
    multiplier = 1.0E0_C_DOUBLE
write(*,*) "... then left boudns shoud be -ive, right +ive"
  ELSE
    multiplier = -1.0E0_C_DOUBLE
write(*,*) "... then left boudns shoud be +ive, right -ive"
  END IF  
  ! FIND the function value of the starting point (SP)
  ! NOTE: findImkM must be called with the parameters from the COMMON block.
  CALL findImkM(i, startZero, SPvalue, df, m)
write(*,*) "Start, and its fn value:", startZero, SPvalue
      ! The fn value at the starting point, so we know which way to search

  ! ************** LOWER BOUND
  ! If fn value at SP is *positive*, only need to creep to the right
  boundL = startZero
  
  IF (SPvalue .LE. 0.0E0_C_DOUBLE) THEN
write(*,*) "SP FUNCTION VALUE NEGATIVE"
    keepSearching = .TRUE.
    DO WHILE (keepSearching)
      ! If fn value at SP is negative, take bold steps left to find lower bound
      boundL = boundL / 2.0E0_C_DOUBLE
      
    CALL findImkM(i, boundL, valueL, df, m)
      
      IF (valueL .GT. 0.0E0_C_DOUBLE ) THEN
        ! Found a lower bound where the fn value is positive
        keepSearching = .FALSE.
      END IF 
    END DO
  END IF

  ! Now creep to the right from boundL (refine the lower bound)
  keepSearching = .TRUE.
  DO WHILE (keepSearching)
write(*,*) "Keep searching.... creeping"
    oldBoundL = boundL
    boundL = boundL * 1.10E0_C_DOUBLE
    CALL findImkM(i, boundL, valueL, df, m)
write(*,*) "NEW BOUND ANd FN VALUE:", boundL, valueL

    IF (valueL .LT. 0.0E0_C_DOUBLE ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundL = oldBoundL
    END IF
  END DO
  zeroL = boundL

  ! ************** UPPER BOUND
  ! SPvaluecomputed above
  boundR = startZero

  ! If fn value at SP is *negative*, only need to creep to the left
  IF (SPvalue .GT. 0.0E0_C_DOUBLE) THEN
    boundR = startZero
    keepSearching = .TRUE.
    
    DO WHILE (keepSearching)
      ! If fn value at SP is positive, take bold steps right to find upper bound
      boundR = boundR * 1.5E0_C_DOUBLE

      CALL findImkM(i, boundR, valueR, df, m)

      IF (valueR .LT. 0.0E0_C_DOUBLE ) THEN
        ! Found an upper bound where the fn value is negative
        keepSearching = .FALSE.
      END IF 
    END DO
  END IF

  ! Now creep to the left from boundR (refine the upper bound)
  keepSearching = .TRUE.
  DO WHILE (keepSearching)
    oldBoundR = boundR
    boundR = boundR * 0.90E0_C_DOUBLE
    CALL findImkM(i, boundR, valueR, df, m)
    IF (valueR .GT. 0.0E0_C_DOUBLE ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundR = oldBoundR
    END IF
  END DO
  zeroR = boundR

END SUBROUTINE improveKZeroBounds
