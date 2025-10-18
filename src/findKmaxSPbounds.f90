SUBROUTINE findKmaxSPbounds(i, startTKmax, kmaxL, kmaxR)
  USE tweedie_params_mod, ONLY: Cphi, Cmu, Cy
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! Inputs
  REAL(KIND=C_DOUBLE), INTENT(IN)    :: startTKmax
  INTEGER(C_INT), INTENT(IN)         :: i
  ! Outputs
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: kmaxL, kmaxR

  ! --- Local Variables ---
  REAL(KIND=C_DOUBLE)    :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE)    :: boundL, boundR, slope, SPslope
  REAL(KIND=C_DOUBLE)    :: oldBoundL, oldBoundR
  LOGICAL         :: keepSearching
  
  EXTERNAL findImkd
  

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! FIND the slope of the starting point (SP)
  ! NOTE: findImkd must be called with the parameters from the COMMON block.
  CALL findImkd(i, startTKmax, SPslope)

  ! ************** LOWER BOUND
  ! If slope at SP is *positive*, only need to creep to the right
  boundL = startTKmax
  
  IF (SPslope <= 0.0d0) THEN
    keepSearching = .TRUE.
    DO WHILE (keepSearching)
      ! If slope at SP is negative, take bold steps left to find lower bound
      boundL = boundL / 2.0d0 
      
      CALL findImkd(i, boundL, slope)
      
      IF (slope > 0.0d0 ) THEN
        ! Found a lower bound where the slope is positive
        keepSearching = .FALSE.
      END IF 
    END DO
  END IF

  ! Now creep to the right from boundL (refine the lower bound)
  keepSearching = .TRUE.
  DO WHILE (keepSearching)
    oldBoundL = boundL
    boundL = boundL * 1.10d0
    CALL findImkd(i, boundL, slope)

    IF (slope < 0.0d0 ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundL = oldBoundL
    END IF
  END DO
  kmaxL = boundL

  ! ************** UPPER BOUND
  ! Re-evaluate SPslope (redundant but safe)
  CALL findImkd(i, startTKmax, SPslope)
  boundR = startTKmax

  ! If slope at SP is *negative*, only need to creep to the left
  IF (SPslope > 0.0d0) THEN
    boundR = startTKmax
    keepSearching = .TRUE.
    
    DO WHILE (keepSearching)
      ! If slope at SP is positive, take bold steps right to find upper bound
      boundR = boundR * 1.5d0

      CALL findImkd(i, boundR, slope)

      IF (slope < 0.0d0 ) THEN
        ! Found an upper bound where the slope is negative
        keepSearching = .FALSE.
      END IF 
    END DO
  END IF

  ! Now creep to the left from boundR (refine the upper bound)
  keepSearching = .TRUE.
  DO WHILE (keepSearching)
    oldBoundR = boundR
    boundR = boundR * 0.90d0
    CALL findImkd(i, boundR, slope)

    IF (slope > 0.0d0 ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundR = oldBoundR
    END IF
  END DO
  kmaxR = boundR
  ! Removed write(*,*)
  
END SUBROUTINE findKmaxSPbounds
