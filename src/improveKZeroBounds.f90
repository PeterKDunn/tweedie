SUBROUTINE improveKZeroBounds(i, m, startTKmax, kmaxL, kmaxR)
  USE tweedie_params_mod, ONLY: Cphi, Cmu, Cy
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Inputs
  REAL(KIND=C_DOUBLE), INTENT(IN)    :: startTKmax
  INTEGER(C_INT), INTENT(IN)         :: i, m
  ! Outputs
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: kmaxL, kmaxR

  ! --- Local Variables ---
  REAL(KIND=C_DOUBLE)    :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE)    :: boundL, boundR, slopeL, slopeR, SPslope
  REAL(KIND=C_DOUBLE)    :: oldBoundL, oldBoundR, f
  LOGICAL         :: keepSearching
  
  EXTERNAL findImkd
  

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! FIND the slope of the starting point (SP)
  ! NOTE: findImkd must be called with the parameters from the COMMON block.
    CALL findImkM(i, boundL, f, slopeL, m)

  ! ************** LOWER BOUND
  ! If slope at SP is *positive*, only need to creep to the right
  boundL = startTKmax
  
  IF (SPslope .LE. 0.0E0_C_DOUBLE) THEN
    keepSearching = .TRUE.
    DO WHILE (keepSearching)
      ! If slope at SP is negative, take bold steps left to find lower bound
      boundL = boundL / 2.0E0_C_DOUBLE
      
    CALL findImkM(i, boundL, f, slopeL, m)
      
      IF (slopeL .GT. 0.0E0_C_DOUBLE ) THEN
        ! Found a lower bound where the slope is positive
        keepSearching = .FALSE.
      END IF 
    END DO
  END IF

  ! Now creep to the right from boundL (refine the lower bound)
  keepSearching = .TRUE.
  DO WHILE (keepSearching)
    oldBoundL = boundL
    boundL = boundL * 1.10E0_C_DOUBLE
    CALL findImkM(i, boundL, f, slopeL, m)

    IF (slopeL .LT. 0.0E0_C_DOUBLE ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundL = oldBoundL
    END IF
  END DO
  kmaxL = boundL

  ! ************** UPPER BOUND
  ! Re-evaluate SPslope (redundant but safe)
  CALL findImkM(i, boundR, f, slopeR, m)
  boundR = startTKmax

  ! If slope at SP is *negative*, only need to creep to the left
  IF (SPslope .GT. 0.0E0_C_DOUBLE) THEN
    boundR = startTKmax
    keepSearching = .TRUE.
    
    DO WHILE (keepSearching)
      ! If slope at SP is positive, take bold steps right to find upper bound
      boundR = boundR * 1.5E0_C_DOUBLE

      CALL findImkM(i, boundR, f, slopeR, m)

      IF (slopeR .LT. 0.0E0_C_DOUBLE ) THEN
        ! Found an upper bound where the slope is negative
        keepSearching = .FALSE.
      END IF 
    END DO
  END IF

  ! Now creep to the left from boundR (refine the upper bound)
  keepSearching = .TRUE.
  DO WHILE (keepSearching)
    oldBoundR = boundR
    boundR = boundR * 0.90E0_C_DOUBLE
    CALL findImkM(i, boundR, f, slopeR, m)
    IF (slopeR .GT. 0.0E0_C_DOUBLE ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundR = oldBoundR
    END IF
  END DO
  kmaxR = boundR

END SUBROUTINE improveKZeroBounds
