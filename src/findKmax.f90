SUBROUTINE findKmax(i, kmax, tmax, mmax, mfirst, leftOfMax)
  ! Finds the value of Kmax, Tmax, and Mmax.
  ! Also return the first value of m (mfirst) amd whether this is to the left of the max (leftOfMax).

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: kmax, tmax
  INTEGER(C_INT), INTENT(OUT)         :: mmax, mfirst, leftOfMax
  INTEGER(C_INT), INTENT(IN)          :: i

  REAL(KIND=C_DOUBLE)     :: pi, t_Start_Point, slope_At_Zero, Imk_value
  REAL(KIND=C_DOUBLE)     :: aimrerr, tmaxL, tmaxR, omega_SP, ratio
  REAL(KIND=C_DOUBLE)     :: current_y, current_mu, current_phi
  LOGICAL(C_BOOL)         :: root_Found
  

  INTERFACE

    FUNCTION findKmaxSP(j) 
      ! Template tehfunction for finding a starting point for finding Kmax
      
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      
      IMPLICIT NONE  
      
      REAL(KIND=C_DOUBLE)   :: findKmaxSP
      INTEGER, INTENT(IN)   :: j
    END FUNCTION findKmaxSP

    
    SUBROUTINE funcd_signature(i_in, t, f, df) BIND(C)
      ! Template for the function for which zeros are sought
      
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      
      IMPLICIT NONE
      
      INTEGER(C_INT), INTENT(IN) :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN) :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
    END SUBROUTINE funcd_signature


    SUBROUTINE rtsafe(i_in, funcd, x1, x2, xacc, root, root_Found) 
      ! Find zeros using (moodified) Newton's method with bisection
      
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
      PROCEDURE(funcd_signature) :: funcd 
    
      INTEGER(C_INT), INTENT(IN)        :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: x1, x2, xacc
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root
      LOGICAL(C_BOOL), INTENT(OUT)      :: root_Found
    END SUBROUTINE rtsafe
      

    SUBROUTINE rtnewton(i_in, funcd, xstart, xacc, root, root_Found) 
      ! Find zeros using (moodified) Newton's method
      
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
      PROCEDURE(funcd_signature)        :: funcd
      
      INTEGER(C_INT), INTENT(IN)        :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: xstart, xacc
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root
      LOGICAL(C_BOOL), INTENT(OUT)      :: root_Found
    END SUBROUTINE rtnewton
      

    SUBROUTINE evaluateImkdZero(i_in, t, f, df) 
      ! Evaluate Im k'(t)  and  Im k''(t)  for solving for Kmax (i.e., Im k'(t) = 0)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      
      IMPLICIT NONE

      INTEGER(C_INT), INTENT(IN)        :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    END SUBROUTINE evaluateImkdZero
      

    SUBROUTINE evaluateImk(i_in, t_in, kmax_out) 
      ! Find Im k(t)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE

      INTEGER(C_INT), INTENT(IN) :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN) :: t_in
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: kmax_out
    END SUBROUTINE evaluateImk
      
  END INTERFACE


  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! --- Initialization ---
  aimrerr = 1.0E-09_C_DOUBLE
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

  ! Find starting points
  CALL evaluateImkd(i, 0.0_C_DOUBLE, slope_At_Zero)

  IF (slope_At_Zero .LE. 0.0_C_DOUBLE) THEN
    ! Im k(t) initially heads downwards
    ! This includes the case 'IF y >= mu'
    ! Nothing to do; easy-peasy:
    mmax = 0
    mfirst = -1
    kmax = 0.0_C_DOUBLE
    tmax = 0.0_C_DOUBLE
    leftOfMax = 0
  ELSE
    ! CASE: IF slope is initially UPWARDS: trickier, esp. with 1 < p < 2
    ! Good starting point often needed

    ratio = current_y / current_mu
  
    omega_SP = -1
    t_Start_Point = current_mu ** (1.0_C_DOUBLE - Cp) * DTAN(omega_SP) /   &
                    ( ( 1.0_C_DOUBLE - Cp) * current_phi)
    IF (ratio .LT. 0.1_C_DOUBLE) THEN
      t_Start_Point = current_mu ** (1.0_C_DOUBLE - Cp) / current_phi *    &
                      DSQRT( 2 * (1.0_C_DOUBLE - ratio))
    END IF
    IF (ratio .GT. 0.9_C_DOUBLE) THEN
      omega_SP = -0.01_C_DOUBLE
      t_Start_Point = current_mu ** (1.0_C_DOUBLE - Cp) * DTAN(omega_SP) /   &
                      ( ( 1.0_C_DOUBLE - Cp) * current_phi)
    END IF

    ! Now find kmax and tmax
    IF (Cpsmall) THEN
      tmaxL = 0.0_C_DOUBLE       ! Since Left bound can be zero
      tmaxR = t_Start_Point * 2.0_C_DOUBLE

      CALL improveKmaxSPBounds(i, t_Start_Point, tmaxL, tmaxR)
      ! Crudely improve the bounds that bracket the starting point for finding Kmax.
      CALL rtsafe(i,                &
                  evaluateImkdZero,   &
                  tmaxL,          &
                  tmaxR,          &
                  aimrerr,        &
                  tmax,           &
                  root_Found)
    ELSE
      ! p > 2
      t_Start_Point = current_mu ** (1.0_C_DOUBLE - Cp) * DTAN(omega_SP) /   &
                      ( ( 1.0_C_DOUBLE - Cp) * current_phi)
      CALL rtnewton(i,              &
                    evaluateImkdZero,   &
                    t_Start_Point,  &
                    aimrerr,        &
                    tmax,           &
                    root_Found)
    END IF

    ! Find mmax, which depends on whether we are working with the PDF or the CDF.
    ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
    ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m pi/y;
    !       the CDF has integrand zeros at Im k(t) =        m pi/y.
    CALL evaluateImk(i, tmax, kmax)
    
    IF (Cpdf) THEN
      mmax = FLOOR(2.0_C_DOUBLE * kmax / pi)
    ELSE
      mmax = FLOOR(kmax / pi)
    END IF

    ! Establish the first value of m to use, and whether the first zero is to the left of kmax
    IF (mmax .GT. 0) THEN
      mfirst = 1
      leftOfMax = 1
    ELSE
      IF (mmax .EQ. 0 ) THEN
        mfirst = 0
        leftOfMax = 0
      ELSE
        ! That is, mmax is LESS THAN 0
        mfirst = -1
        leftOfMax = 0
      ENDIF 
    END IF
  END IF

  CONTAINS
    
    SUBROUTINE improveKmaxSPBounds(i, startTKmax, tmaxL, tmaxR)
      ! Crudely improve the bounds that bracket the starting point for finding Kmax.
      ! Sometime a decent starting point is crucial for timely convergence.
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
    
      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: startTKmax
      INTEGER(C_INT), INTENT(IN)          :: i
      REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: tmaxL, tmaxR
    
      ! --- Local Variables ---
      REAL(KIND=C_DOUBLE)     :: current_y, current_mu, current_phi
      REAL(KIND=C_DOUBLE)     :: boundL, boundR, slopeL, slopeR
      REAL(KIND=C_DOUBLE)     :: oldBoundL, oldBoundR
      INTEGER(C_INT)          :: max_Search, search_Its
      LOGICAL                 :: keep_Searching
      
      EXTERNAL evaluateImkd
      
    
      ! Grab the relevant scalar values for this iteration:
      current_y    = Cy(i)    ! Access y value for index i
      current_mu   = Cmu(i)   ! Access mu value for index i
      current_phi  = Cphi(i)  ! Access phi value for index i
    
      max_Search = 10
      
      !!!!! LOWER BOUND
      !   - If slope at SP is *positive* (which it should be), only need to creep to the right
      boundL = tmaxL
      boundR = tmaxR

      CALL evaluateImkd(i, boundL, slopeL)
      CALL evaluateImk(i, boundL, Imk_value)
      
      ! Slope at starting point should always be positive if we are here:
      ! if the slope is negative, Im k(t) heads down and m = -1, -2, ...

      ! - Now creep to the right from boundL (refine the lower bound)
      keep_Searching = .TRUE.
      DO WHILE (keep_Searching)
        oldBoundL = boundL
        boundL = (boundL + 1.0E-2_C_DOUBLE) * 1.250E0_C_DOUBLE
        CALL evaluateImkd(i, boundL, slopeL)
        CALL evaluateImk(i, boundL, Imk_value)

        IF ( (slopeL .LT. 0.0E0_C_DOUBLE ) .OR.   &
             (Imk_value .LT. 0.0_C_DOUBLE) ) THEN
          ! - Gone too far, so keep previous bound
          keep_Searching = .FALSE.
          boundL = oldBoundL
        END IF
      END DO
      tmaxL = boundL


      !!!!! UPPER BOUND
      
      ! If the improved lower bouond is LARGER than the original upper bound... FIX! 
      IF (tmaxL .GT. boundR) boundR = tmaxL * 2.0_C_DOUBLE 
      
      ! Find the slope at boundR the starting point (SP)
      CALL evaluateImkd(i, boundR, slopeR)
      
      ! Find the value if kmax at the starting point (SP)
      CALL evaluateImk(i, boundR, Imk_value)

      ! A valid upper bound:
      ! - must have a negative SLOPE (must be negative)i.e., heading down), *AND*
      ! - must have a positive vale of kmax (or we may have a local maximum only) 
    
      search_Its = 0
      ! - If slope at SP is *positive*, and kmax positive, need bold steps to the right
      IF ( (slopeR .GT. 0.0E0_C_DOUBLE) .AND.  &
           (Imk_value .GT. 0.0_C_DOUBLE) ) THEN
        boundR = startTKmax
        keep_Searching = .TRUE.
        
        DO WHILE (keep_Searching)
          ! - If slope at SP is positive, take bold steps right to find upper bound
          search_Its = search_Its + 1
          boundR = (boundR + 0.1_C_DOUBLE) * 2.0E0_C_DOUBLE
    
          CALL evaluateImkd(i, boundR, slopeR)
          CALL evaluateImk(i, boundR, Imk_value)

          IF ( (slopeR .LT. 0.0E0_C_DOUBLE ) .AND.  &
               (Imk_value .GT. 0.0_C_DOUBLE) ) THEN
            ! - Found an upper bound where the slope is negative, and kmax is positive
            keep_Searching = .FALSE.
          END IF 
          IF ( search_Its .GT. max_Search) keep_Searching = .FALSE.
        END DO
      END IF

      ! - Now creep to the left from boundR (refine the upper bound)
      keep_Searching = .TRUE.
      DO WHILE (keep_Searching)
        oldBoundR = boundR
        boundR = boundR * 0.90E0_C_DOUBLE
        CALL evaluateImkd(i, boundR, slopeR)
        CALL evaluateImk(i, boundR, Imk_value)

        IF ( (slopeR .GT. 0.0E0_C_DOUBLE) .AND.   &
             (Imk_value .GT. 0.0_C_DOUBLE) ) THEN
          ! - Gone too far, so keep previous bound
          keep_Searching = .FALSE.
          boundR = oldBoundR
        END IF
      END DO
      tmaxR = boundR

    END SUBROUTINE improveKmaxSPBounds


END SUBROUTINE findKmax
