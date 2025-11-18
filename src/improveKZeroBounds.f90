SUBROUTINE improveKZeroBounds(i, m, leftOfMax, zeroMid, zeroL, zeroR)
  ! Improve the bounds that bracket the zero of Im k(t).
  ! A decent starting point is sometimes crucial to timely convergence.
  
  USE tweedie_params_mod, ONLY: Cphi, Cmu, Cy
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE Calcs_Imag, ONLY: evaluateImkM
  
  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroMid
  INTEGER(C_INT), INTENT(IN)          :: i, m
  REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroL, zeroR
  LOGICAL(C_BOOL), INTENT(IN)         :: leftOfMax

  REAL(KIND=C_DOUBLE)     :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE)     :: valueL, valueR, multiplier
  REAL(KIND=C_DOUBLE)     :: df, valueMid
  INTEGER(C_INT)          :: maxSearch


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i



  ! Initialisation
  maxSearch = 10  ! Don't spend too long, so set limit

  ! Set multipier: this adjust the sign depending on whether we are
  ! left of the max (so left bound is negative) or to the right of
  ! the max (so left bound is positive)
  IF (leftOfMax) THEN
    multiplier = -1.0E0_C_DOUBLE
  ELSE
    multiplier = 1.0E0_C_DOUBLE
  END IF


  ! TWO ROLES:
  ! a) if the bounds actually do bound the zero (i.e, bound give opposite signs for the function
  !    for which zeros are sought): improve if we can, Smetimes, a good start pt is necessary.
  ! b) if the bounds actually do NOT bound the zero: find bounds that do!

  ! FIND the function value of the starting point (SP)
  zeroMid = (zeroL + zeroR) / 2.0_C_DOUBLE
  CALL evaluateImkM(i, zeroMid, valueMid, df, m)
  CALL evaluateImkM(i, zeroL, valueL, df, m)
  CALL evaluateImkM(i, zeroR, valueR, df, m)


  ! CHECK IF BOUNDS REALLY DO BOUND THE ZERO:
  IF ( (valueL * valueR) .GE. 0.0_C_DOUBLE) THEN
    ! Bounds DO NOT trap the zero

    ! The solution depends on what side of the max we are.
    ! If to the LEFT of the max of Im k(t), zeroL should give a -ive value; zeroR a +ive value.
    IF ( leftOfMax) THEN
      DO WHILE ( valueL .GT. 0.0_C_DOUBLE) 
        ! We are on the LEFT of the maximum of Im k(t), but the L bound gives a +ive value.
        ! So we need to go LEFT a little.
        zeroL = (zeroL - 0.1_C_DOUBLE) * 0.95_C_DOUBLE
        CALL evaluateImkM(i, zeroL, valueL, df, m)
      END DO



      ! The solution depends on what side of the max we are.
      ! If to the LEFT of the max of Im k(t), zeroL should give a -ive value; zeroR a +ive value.
      DO WHILE ( valueR .LT. 0.0_C_DOUBLE) 
        ! We are on the LEFT of the maximum of Im k(t), but the R bound gives a -ive value.
        ! So we need to go RIGHT a little.
        zeroR = (zeroR + 0.1_C_DOUBLE) * 1.05_C_DOUBLE
        CALL evaluateImkM(i, zeroR, valueR, df, m)
      END DO
    END IF





    ! The solution depends on what side of the max we are.
    ! If to the RIGHT of the max of Im k(t), zeroL should give a +ive value; zeroR a -ive value.
    IF ( .NOT.(leftOfMax) ) THEN
      DO WHILE ( valueL .LT. 0.0_C_DOUBLE) 
        ! We are on the RIGHT of the maximum of Im k(t), but the L bound gives a -ive value.
        ! So we need to go LEFT a little.
        zeroL = (zeroL - 0.1_C_DOUBLE) * 0.95_C_DOUBLE
        CALL evaluateImkM(i, zeroL, valueL, df, m)
      END DO



      ! The solution depends on what side of the max we are.
      ! If to the RIGHT of the max of Im k(t), zeroL should give a -ive value; zeroR a +ive value.
      DO WHILE ( valueR .GT. 0.0_C_DOUBLE) 
        ! We are on the RIGHT of the maximum of Im k(t), but the R bound gives a +ive value.
        ! So we need to go RIGHT a little.
        zeroR = (zeroR + 0.1_C_DOUBLE) * 1.05_C_DOUBLE
        CALL evaluateImkM(i, zeroR, valueR, df, m)
      END DO
    END IF




  ELSE
    ! Bounds DO trap the zero, so improve a little
  
    ! Find a point halfway between bounds.
    ! If the new point has same sign as L/R bound, make that the new L/R bounds.
    IF ( (valueMid * valueL) .GE. 0.0_C_DOUBLE) THEN
      zeroL = zeroMid
    ELSE IF ( (valueMid * valueR) .GE. 0.0_C_DOUBLE) THEN
      zeroR = zeroMid
    END IF
      
    ! And once more ONLY
    zeroMid = (zeroL + zeroR) / 2.0_C_DOUBLE
    CALL evaluateImkM(i, zeroMid, valueMid, df, m)
    CALL evaluateImkM(i, zeroL, valueL, df, m)
    CALL evaluateImkM(i, zeroR, valueR, df, m)

    ! Find a point halfway between bounds.
    ! If the new point has same sign as L/R bound, make that the new L/R bounds.
    IF ( (valueMid * valueL) .GE. 0.0_C_DOUBLE) THEN
      zeroL = zeroMid
    ELSE IF ( (valueMid * valueR) .GE. 0.0_C_DOUBLE) THEN
      zeroR = zeroMid
    END IF
    
    zeroMid = (zeroL + zeroR) / 2.0_C_DOUBLE

  END IF
  
  RETURN

END SUBROUTINE improveKZeroBounds
