SUBROUTINE findKmax(i, kmax, tmax, mmax, mfirst, leftOfMax)
  ! Finds the value of Kmax, Tmax, and Mmax.
  ! Also return the first value of m (mfirst) amd whether this is to the left of the max (leftOfMax).

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE Calcs_Imag
  USE Calcs_Solvers
  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: kmax, tmax
  INTEGER(C_INT), INTENT(OUT)         :: mmax, mfirst
  LOGICAL(C_BOOL), INTENT(OUT)        :: leftOfMax
  INTEGER(C_INT), INTENT(IN)          :: i

  REAL(KIND=C_DOUBLE)     :: pi, t_Start_Point, slope_At_Zero, Imk_value
  REAL(KIND=C_DOUBLE)     :: aimrerr, tmaxL, tmaxR, ratio, threshold, t_small
  REAL(KIND=C_DOUBLE)     :: current_y, current_mu, current_phi
  LOGICAL(C_BOOL)         :: error
  

  ! A. If Im k(t) heads down initially: easy: kmax = tmax = mmax = 0
  ! B. Otherwise, we need to find kmax etc by solving Im k'(t) = 0
  !    with a suitable starting point:
  !    1. If tmax appears that it will be very large, find approx start pt
  !    2. Otherwise find start point numerically:
  !       a. If p near 1, approximate starting point by p=1 case analytically
  !       b. Otherwise, use standard starting point


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! --- Initialization ---
  aimrerr = 1.0E-09_C_DOUBLE
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  threshold = 1.0E5_C_DOUBLE ! Threshold for a "large" tmax


  ! Find starting points
  CALL evaluateImkd(i, 0.0_C_DOUBLE, slope_At_Zero)
!WRITE(*,*) "INITIAL SLOPE:", slope_At_Zero

  IF (slope_At_Zero .LE. 0.0_C_DOUBLE) THEN
!WRITE(*,*) "SLOPE LE 0"
    ! A. Im k(t) initially heads DOWNWARDS
    !
    !    This includes the case 'IF y >= mu'
    !    Nothing to do; easy-peasy:
    mmax = 0
    mfirst = -1
    kmax = 0.0_C_DOUBLE
    tmax = 0.0_C_DOUBLE
    leftOfMax = .FALSE.
    
    RETURN
  ELSE
!WRITE(*,*) "SLOPE > 0"
    ! B. CASE: IF slope is initially UPWARDS: trickier, esp. with 1 < p < 2
    !    Need to solve  Im k'(t) = 0.

    IF ( ( current_y ** (1.0_C_DOUBLE - Cp) / current_phi) .GT.  &
         (10.0_C_DOUBLE * threshold ) ) THEN 
!WRITE(*,*) "APPROX SP: very very Large t"
      ! B.1. Sometimes kmax and co are MASSIVE, making things slow and difficult, and 
      !      sometimes convergence is difficult also.
      ! In some cases, tmax, kmax and mmax are HUGE; e.g., when we have y=0.001, 
      !     power=6, mu=.01, phi=.01, verbose=TRUE), we can see that 
      ! tmax > 1,000,000, kmax > 80,000 and mmax > 25,000...
      !
      ! When p near 1, we can solve Im k'(t) = 0 when p = 1:
      !    t = acos(y/mu)/phi. 
      ! Much more accurate!
      !
      ! When p near 2, do similar with p = 2; need to check, but I get:
      !    t = 1/(phi*mu) * sqrt( (mu - y)/y )
      ! provided mu - y (which checks out).
      !
      ! This hopefully will flag potentially problematically large kmax/tmax/mmax:

      tmax = threshold
      CALL evaluateImk(i, t_Start_Point, kmax, error)    ! Now find the corresponding value of kmax
      IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, t_Start_Point, 1)

    ELSE
      ! Compute ratio
      ratio = current_y / current_mu
      
      ! 1. Try small-t approximation first (for ratio < 1)
      IF (ratio < 1.0_C_DOUBLE) THEN
          t_small = (current_mu**(1.0_C_DOUBLE - Cp) / current_phi) *    &
                    DSQRT( (2.0_C_DOUBLE/Cp) * (1.0_C_DOUBLE - ratio) )
          ! Clamp tiny values
          IF (t_small < 1.0E-12_C_DOUBLE) t_small = 1.0E-12_C_DOUBLE
          t_Start_Point = t_small
      ELSE
          ! 2. Small-t not applicable → standard start
          IF (Cp < 1.1) THEN
              t_Start_Point = DACOS(current_y/current_mu)/current_phi
          ELSE
              ! Use large-t blended guess
              t_Start_Point = tInitialGuess()   ! This computes large-t approx
              IF (t_Start_Point <= 0.0_C_DOUBLE) t_Start_Point = 1.0E-12_C_DOUBLE
          END IF
      END IF
!WRITE(*,*) "GOT START PT:", t_Start_Point
      ! Now find kmax and tmax using this t_Start_Point
!      IF (Cpsmall) THEN
        tmaxL = 0.0_C_DOUBLE       ! Since Left bound can be zero
        tmaxR = t_Start_Point * 2.0_C_DOUBLE
!WRITE(*,*) "  BOUNDS 1: RTSAFE", tmaxL, tmaxR
        CALL improveKmaxSPBounds(i, t_Start_Point, tmaxL, tmaxR)
!WRITE(*,*) "  BOUNDS 1 revises: RTSAFE", tmaxL, tmaxR
        ! Crudely improve the bounds that bracket the starting point for finding Kmax.
        CALL rtsafe(i,                  &
                      evaluateImkdZero,   &
                      tmaxL,              &
                      tmaxR,              &
                      aimrerr,            &
                      tmax,               &
                      error)
        IF (error) CALL DBLEPR("ERROR: cannot solve", -1, tmax, 1)
    END IF  

    ! Find mmax, which depends on whether we are working with the PDF or the CDF.
    ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
    ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m pi/y;
    !       the CDF has integrand zeros at Im k(t) =        m pi/y.

    ! Find kmax from tmax
    CALL evaluateImk(i, tmax, kmax, error)
    IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, tmax, 1)

    ! Find mmax from kmax
    IF (Cpdf) THEN
      mmax = FLOOR(2.0_C_DOUBLE * kmax / pi)
    ELSE
      mmax = FLOOR(kmax / pi)
    END IF

    ! Establish the first value of m to use, and whether the first zero is to the left of kmax
    IF (mmax .GT. 0) THEN
      mfirst = 1
      leftOfMax = .TRUE.
    ELSE
      IF (mmax .EQ. 0 ) THEN
        mfirst = 0
        leftOfMax = .FALSE.
      ELSE
        ! That is, mmax is LESS THAN 0
        mfirst = -1
        leftOfMax = .FALSE.
      ENDIF 
    END IF
  END IF


  CONTAINS
    
    
    FUNCTION tInitialGuess() RESULT(t0)
    
      IMPLICIT NONE
  
      REAL(C_DOUBLE)              :: t0, pi
      REAL(C_DOUBLE)              :: ratio, r, t_small, t_large, c, angle
    
      pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    
      ! Weighting exponent for interpolation
      r = 4.0E0_C_DOUBLE ! Artibtrary
    
      ! Ratio y/mu
      ratio = current_y / current_mu
    
      ! --- Small-t approximation ---
      IF (ratio < 1.0E0_C_DOUBLE) THEN
         t_small = (current_mu**(1.0E0_C_DOUBLE - Cp) / current_phi) *      &
                    DSQRT( (2.0E0_C_DOUBLE/Cp) * (1.0E0_C_DOUBLE - ratio) )
!WRITE(*,*) "ratio=", ratio
!WRITE(*,*) "entering small-t branch"

      ELSE
         ! If y >= mu, no small-t root: set small approx to zero
         t_small = 0.0E0_C_DOUBLE
      ENDIF
    
      ! --- Large-t approximation ---
      angle   = 0.5D0 * pi / (Cp - 1.0E0_C_DOUBLE)
      c       = DCOS(angle)**(Cp - 1.0E0_C_DOUBLE)
     t_large = ( current_mu**(2.0_C_DOUBLE*(1.0_C_DOUBLE - Cp)) /   &
                ((1.0_C_DOUBLE - Cp) * current_phi) ) *             &
                ( current_y**(Cp - 1.0_C_DOUBLE) )
    
      ! --- Blended initial guess ---
      t0 = ratio**r * t_small + (1.0E0_C_DOUBLE - ratio**r) * t_large
    
    END FUNCTION tInitialGuess

    
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
      LOGICAL(C_BOOL)         :: keep_Searching, error
      
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
      CALL evaluateImk(i, boundL, Imk_value, error)
      IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, boundL, 1)
      
      ! Slope at starting point should always be positive if we are here:
      ! if the slope is negative, Im k(t) heads down and m = -1, -2, ...

      ! - Now creep to the right from boundL (refine the lower bound)
      keep_Searching = .TRUE.
      DO WHILE (keep_Searching)
        oldBoundL = boundL
        boundL = (boundL + 1.0E-2_C_DOUBLE) * 1.250E0_C_DOUBLE
        CALL evaluateImkd(i, boundL, slopeL)
        CALL evaluateImk(i, boundL, Imk_value, error)
        IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, boundL, 1)

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
      CALL evaluateImk(i, boundR, Imk_value, error)
      IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)

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
          CALL evaluateImk(i, boundR, Imk_value, error)
          IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)

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
        CALL evaluateImk(i, boundR, Imk_value, error)
        IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)

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
