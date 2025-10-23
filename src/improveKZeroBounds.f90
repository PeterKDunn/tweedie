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
  INTEGER(C_INT)          :: maxSearch, itsSearch
  LOGICAL                 :: keepSearching
  
  EXTERNAL findImkM
  

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  maxSearch = 10
  ! Set multipier: this adjust the sign depending on whether we are
  ! left of the max (so left bound is negative) or to the right of
  ! the max (so left bound is positive)
  IF (leftOfMax .EQ. 1) THEN
    multiplier = -1.0E0_C_DOUBLE
  ELSE
    multiplier = 1.0E0_C_DOUBLE
  END IF

  ! FIND the function value of the starting point (SP)
  ! NOTE: findImkM must be called with the parameters from the COMMON block.
  CALL findImkM(i, startZero, SPvalue, df, m)

      ! The fn value at the starting point, so we know which way to search

  ! ************** LOWER BOUND
  ! If fn value at SP is *positive*, only need to creep to the right
  boundL = startZero
  itsSearch = 0
  
  IF ( (multiplier * SPvalue) .LE. 0.0E0_C_DOUBLE) THEN

    itsSearch = itsSearch + 1
    keepSearching = .TRUE.
    DO WHILE (keepSearching)
      ! If fn value at SP is negative, take bold steps left to find lower bound
      boundL = boundL / 1.50E0_C_DOUBLE
      
    CALL findImkM(i, boundL, valueL, df, m)
      
      IF ( (multiplier * valueL) .GT. 0.0E0_C_DOUBLE ) THEN
        ! Found a lower bound where the fn value is positive
        keepSearching = .FALSE.
      END IF 
    END DO
  END IF

  ! Now creep to the right from boundL (refine the lower bound)
  keepSearching = .TRUE.
  itsSearch = 0
  DO WHILE (keepSearching)
    itsSearch = itsSearch + 1

    oldBoundL = boundL
    boundL = boundL * 1.10E0_C_DOUBLE
    CALL findImkM(i, boundL, valueL, df, m)

    IF ( (multiplier * valueL) .LT. 0.0E0_C_DOUBLE ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundL = oldBoundL
    END IF
    IF (itsSearch .GT. maxSearch) keepSearching = .FALSE.

  END DO
  zeroL = boundL

  ! ************** UPPER BOUND
  ! SPvaluecomputed above
  boundR = startZero

  itsSearch = 0
  ! If fn value at SP is *negative*, only need to creep to the left
  IF ( (multiplier * SPvalue) .GT. 0.0E0_C_DOUBLE) THEN
    itsSearch = itsSearch + 1
    boundR = startZero
    keepSearching = .TRUE.
    
    DO WHILE (keepSearching)
      ! If fn value at SP is positive, take bold steps right to find upper bound
      boundR = boundR * 1.5E0_C_DOUBLE

      CALL findImkM(i, boundR, valueR, df, m)

      IF ( (multiplier * valueR) .LT. 0.0E0_C_DOUBLE ) THEN
        ! Found an upper bound where the fn value is negative
        keepSearching = .FALSE.
      END IF 
      IF (itsSearch .GT. maxSearch) keepSearching = .FALSE.
      END DO
  END IF

  ! Now creep to the left from boundR (refine the upper bound)
  itsSearch = 0
  keepSearching = .TRUE.
  DO WHILE (keepSearching)
    itsSearch = itsSearch + 1
    oldBoundR = boundR
    boundR = boundR * 0.90E0_C_DOUBLE
    CALL findImkM(i, boundR, valueR, df, m)
    IF ( (multiplier * valueR).GT. 0.0E0_C_DOUBLE ) THEN
      ! Gone too far! Keep previous bound
      keepSearching = .FALSE.
      boundR = oldBoundR
    END IF
    IF (itsSearch .GT. maxSearch) keepSearching = .FALSE.
  END DO
  zeroR = boundR

  RETURN

END SUBROUTINE improveKZeroBounds
