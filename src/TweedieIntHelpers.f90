
MODULE TweedieIntHelpers
  
  USE tweedie_params_mod
  USE Calcs_K
  USE Calcs_Imag
  USE Calcs_Real
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE R_interfaces

  IMPLICIT NONE

CONTAINS 
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE checkStopPreAcc(tmax, zeroL, consecSmallCount, &
                             stop_PreAccelerate, converged_Pre, error, &
                             tDone, areaSoFar, absArea, aimrerr)
    ! Determine if it is OK to stop pre-accelerating, and start using acceleration
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN) :: tmax
    REAL(KIND=C_DOUBLE), INTENT(IN) :: zeroL
    LOGICAL(C_BOOL), INTENT(OUT)    :: stop_PreAccelerate, converged_Pre
    LOGICAL(C_BOOL), INTENT(INOUT)  :: error
    INTEGER(C_INT), INTENT(INOUT)   :: consecSmallCount
    REAL(KIND=C_DOUBLE), INTENT(IN) :: tDone      ! Integration so far covers [0, tDone]
    REAL(KIND=C_DOUBLE), INTENT(IN) :: areaSoFar  ! Current estimate of the integral over [0, tDone]
    REAL(KIND=C_DOUBLE), INTENT(IN) :: absArea    ! Sum of |region areas| over [0, tDone]
    REAL(KIND=C_DOUBLE), INTENT(IN) :: aimrerr    ! Target relative error

    ! Local vars
    INTEGER(C_INT)        :: nmax
    REAL(KIND=C_DOUBLE)   :: MM, Rek, Rekd, tstop, Imk, lambda
    REAL(KIND=C_DOUBLE)   :: condEnvelope, logTail, tailTol
    REAL(KIND=C_DOUBLE)   :: wRe, wIm, eRe, eIm
    LOGICAL(C_BOOL)       :: errorHere
    
    ! NOTE: 
    ! If  stop_PreAccelerate  is  .TRUE.  it means to stop pre-accelerating
    ! If  converged_Pre  is  .TRUE.  it means to convergence has been identified during pre-accelerating
    ! These are NOT necessarily the same.
    !
    ! If convergence is detected (converged_Pre is .TRUE.), then  stop_PreAccelerate  is also .TRUE.
    !
    ! But sometimes, pre-acceleration should stop (stop_PreAccelerate  is  .TRUE.)  even if
    ! convergence is NOT detected; for example:
    !  - the maximum number of iterations has been reached (so stop_PreAccelerate is .TRUE, 
    !    but  converged_Pre  is  .FALSE.); or
    ! - things are so well behaved that we can go straight to acceleration (stop_PreAccelerate is 
    !    .TRUE.  but  converged_Pre  is still  .FALSE.)
    

    ! Initialise
    stop_PreAccelerate = .FALSE.
    converged_Pre = .FALSE.

    ! For p > 2: a rigorous, closed-form upper bound on the remaining tail.
    ! Re k(t) is strictly decreasing and Re k(t) <= C - B t^alpha for all t > 0,
    ! so  |int_{tDone}^inf integrand dt|  <=  int_{tDone}^inf exp(Re k(t)) dt  (PDF),
    ! with an extra factor 1/tDone for the CDF (integrand carries 1/t).
    ! This does not depend on tmax (which, in the 'very large tmax' fallback of
    ! findKmax, is only a placeholder), so it detects convergence long before
    ! zeroL > tmax in cases like p=3, mu=5, phi=0.01, y=0.001.
    ! Because the bound is rigorous and monotone in tDone, one pass is enough:
    ! no consecutive-region guard is needed.
    IF ( (Cp .GT. 2.0_C_DOUBLE) .AND. (tDone .GT. 0.0_C_DOUBLE) ) THEN
      logTail = logTailBoundPgt2(tDone)
      IF (.NOT. Cpdf) logTail = logTail - DLOG(tDone)
      ! Tolerance: relative to the quantity actually returned (on the scale of
      ! the integral), but never below the roundoff level of the summation
      ! (cancellation can leave the result itself at noise level, e.g. when the
      ! true density is ~ 0).
      !   PDF:  result = area/pi                 -> scale |area|
      !   CDF:  result = 0.5 -/+ area/pi (p > 2)  -> scale pi |0.5 -/+ area/pi|
      IF (Cpdf) THEN
        tailTol = aimrerr * MAX( DABS(areaSoFar), 1.0E-16_C_DOUBLE * absArea )
      ELSE
        IF (Ctail) THEN
          tailTol = PI * DABS(0.5_C_DOUBLE + areaSoFar / PI)
        ELSE
          tailTol = PI * DABS(0.5_C_DOUBLE - areaSoFar / PI)
        END IF
        tailTol = aimrerr * MAX( tailTol, 1.0E-16_C_DOUBLE * (absArea + 0.5_C_DOUBLE * PI) )
      END IF
      IF ( tailTol .GT. 0.0_C_DOUBLE ) THEN
        IF ( logTail .LT. DLOG(tailTol) ) THEN
          stop_PreAccelerate = .TRUE.
          converged_Pre = .TRUE.
          RETURN
        END IF
      END IF
    END IF

    ! Stop condition for pre-acceleration.
    ! Ensure that we have passed the peak of Im k(t), so that acceleration can be used
    
    IF ( Cp .GT. 2.0_C_DOUBLE) THEN
      IF (current_y .GT. current_mu) THEN
        ! Stop pre-accelerating, and go straight to acceleration 
        stop_PreAccelerate = .TRUE.
        converged_Pre = .FALSE.
      ELSE
        IF (zeroL .GT. tmax) THEN
          ! past tmax, we can usually stop pre-accelerating, and go straight to acceleration 
          stop_PreAccelerate = .TRUE.
        END IF
      END IF
    ELSE
      IF (current_y .GT. current_mu) THEN
        ! Stop pre-accelerating, and go straight to acceleration 
        stop_PreAccelerate = .TRUE.
      ELSE
        MM = 1.0_C_DOUBLE / (2.0_C_DOUBLE  * (Cp - 1.0_C_DOUBLE))
        IF ( MM .LE. 1.0_C_DOUBLE ) THEN
          stop_PreAccelerate = .TRUE.
        ELSE
          ! Check when t is larger than the last turning point of Re k(t) 
          ! Or when exp{Re(k)/t} is so small that it makes no difference... but 
          ! care is needed: Re k(t) is not necessarily convex here
           IF ( ABS( MM - FLOOR(MM) ) < 1.0E-09_C_DOUBLE ) then
              nmax = FLOOR(MM) - 1_C_INT
           ELSE
              nmax = FLOOR(MM)
           END IF
           tstop = current_mu**(1.0_C_DOUBLE - Cp) / ((1.0_C_DOUBLE - Cp) * current_phi) *   & 
                   DTAN( DBLE(nmax) * PI * (1.0_C_DOUBLE - Cp) )
           IF (zeroL .GT. tstop) stop_PreAccelerate = .TRUE.
        END IF
      END IF
    END IF

    ! Sometimes this takes forever to flag  stop_PreAccelerate  as .TRUE.
    ! so also check if exp{Re k(t)/t} is very small
    CALL evaluateRek( zeroL, Rek, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: cSPreAcc: Rek not found at", -1, zeroL, 1)
      RETURN
    END IF

    CALL evaluateRekd(zeroL, Rekd, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: cSPreAcc: Rekd not found at", -1, zeroL, 1)
      RETURN
    END IF
    
    IF (CpSmall) THEN
      CALL evaluateLambda(lambda)
      CALL evaluateImk(zeroL, Imk, errorHere)
      IF (errorHere) THEN
        error = .TRUE.
        RETURN
      END IF
      ! Use the AMPLITUDE of the oscillation, not the integrand's value here:
      ! zeroL is (close to) a zero of the integrand, so its value there is ~0
      ! whatever the amplitude, which made convergence look reached after a
      ! few regions (e.g. p = 1.05, phi = 0.001: declared converged at t = 22,
      ! where |exp(k(t))| is still ~0.6; F was then wrong by ~2e-2, unflagged).
      !
      ! With K(t) = k(t) + i t y, the numerator of the integrand is
      !   Im or Re of  (exp(K) - exp(-lambda)) * exp(-i t y),
      ! and K(t) -> -lambda as t -> infinity, so the two parts cancel. Its
      ! amplitude is therefore  exp(-lambda) * |exp(w) - 1|,
      !   w = (Re k + lambda) + i (Im k + t y),
      ! computed without cancellation as
      !   Re(exp(w) - 1) = expm1(a) cos(b) - 2 sin^2(b/2),  Im = exp(a) sin(b).
      wRe = Rek + lambda
      wIm = Imk + zeroL * current_y
      eRe = expm1Local(wRe) * DCOS(wIm) - 2.0_C_DOUBLE * DSIN(0.5_C_DOUBLE * wIm)**2
      eIm = DEXP(wRe) * DSIN(wIm)
      condEnvelope = DEXP(-lambda) * DSQRT(eRe*eRe + eIm*eIm)
    ELSE
      condEnvelope = DEXP(Rek)
    END IF

    IF (zeroL .GT. 0.0_C_DOUBLE) THEN

      IF (CpSmall) THEN
        ! For CpSmall, never declare convergence at the loose 1e-7 level --
        ! only at a threshold tight enough to match aimrerr, so pre-acceleration
        ! doesn't exit early and skip the properly-checked acceleration phase.
        IF ( (condEnvelope/zeroL) .LT. 1.0E-07_C_DOUBLE ) THEN
          consecSmallCount = consecSmallCount + 1_C_INT
        ELSE
          consecSmallCount = 0_C_INT
        END IF

        IF ( (consecSmallCount .GE. 3_C_INT) .AND. &
             ((condEnvelope/zeroL) .LT. 1.0E-13_C_DOUBLE) ) THEN
          stop_PreAccelerate = .TRUE.
          converged_Pre = .TRUE.
        END IF

      ELSE
        ! Guard: require several consecutive regions past tmax satisfying
        ! the envelope/slope test before declaring convergence, not just
        ! one. A single region immediately after crossing tmax can satisfy
        ! the loose 1e-7 test by chance for this (p>2, y<mu) regime,
        ! causing catastrophic truncation of the integral.
        IF ( zeroL .GT. tmax ) THEN
          IF ( ( (condEnvelope/zeroL) .LT. 1.0E-07_C_DOUBLE) .AND. (Rekd .LT. 0.0_C_DOUBLE) ) THEN
            consecSmallCount = consecSmallCount + 1_C_INT
          ELSE
            consecSmallCount = 0_C_INT
          END IF

          IF ( consecSmallCount .GE. 3_C_INT ) THEN
            stop_PreAccelerate = .TRUE.
            converged_Pre = .TRUE.
          END IF

          IF ( (condEnvelope/zeroL) .LT. 1.0E-15_C_DOUBLE .AND. &
               (consecSmallCount .GE. 3_C_INT) ) THEN
            stop_PreAccelerate = .TRUE.
            converged_Pre = .TRUE.
          END IF
        END IF
      END IF

    END IF
    
    ! If converged, then always stop pre-accelerating
    IF (converged_Pre) stop_PreAccelerate = .TRUE.

  END SUBROUTINE checkStopPreAcc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  PURE FUNCTION expm1Local(x) RESULT(v)
    ! exp(x) - 1, accurate also for small |x|
    REAL(KIND=C_DOUBLE), INTENT(IN) :: x
    REAL(KIND=C_DOUBLE)             :: v, u
    IF (DABS(x) .GT. 1.0E-5_C_DOUBLE) THEN
      v = DEXP(x) - 1.0_C_DOUBLE
    ELSE
      ! Taylor series: x + x^2/2 + x^3/6 is exact to double precision here
      u = x
      v = u * (1.0_C_DOUBLE + u * (0.5_C_DOUBLE + u / 6.0_C_DOUBLE))
    END IF
  END FUNCTION expm1Local


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  FUNCTION logTailBoundPgt2(t0) RESULT(lb)
    ! log of a rigorous upper bound on  int_{t0}^inf exp(Re k(t)) dt,  valid for p > 2.
    !
    ! With  alpha = (p-2)/(p-1)  in (0,1),  C = mu^(2-p) / (phi (p-2)),
    !       c = (p-1) phi mu^(p-1):
    !   Re k(t) = C - C Re[(1 - i c t)^alpha]  <=  C - B t^alpha,
    !   B = C cos(alpha pi/2) c^alpha,
    ! because |arg(1 - i c t)| < pi/2  and  |1 - i c t| >= c t.
    ! Substituting  x = B t^alpha:
    !   int_{t0}^inf exp(C - B t^alpha) dt = exp(C) Gamma(s, x0) / (alpha B^s),  s = 1/alpha,
    ! and  Gamma(s, x) <= x^(s-1) exp(-x) * x / (x - s + 1)  for  x > s - 1.
    ! Computed on the log scale (B^s can overflow when p is near 2).
    ! Returns HUGE when the bound is not yet usable (x <= s).

    IMPLICIT NONE
    REAL(KIND=C_DOUBLE), INTENT(IN) :: t0
    REAL(KIND=C_DOUBLE)             :: lb
    REAL(KIND=C_DOUBLE)             :: alpha, C, cc, B, s, x

    lb = HUGE(1.0_C_DOUBLE)
    IF ( (Cp .LE. 2.0_C_DOUBLE) .OR. (t0 .LE. 0.0_C_DOUBLE) ) RETURN

    alpha = (Cp - 2.0_C_DOUBLE) / (Cp - 1.0_C_DOUBLE)
    C     = current_mu**(2.0_C_DOUBLE - Cp) / (current_phi * (Cp - 2.0_C_DOUBLE))
    cc    = (Cp - 1.0_C_DOUBLE) * current_phi * current_mu**(Cp - 1.0_C_DOUBLE)
    B     = C * DCOS(0.5_C_DOUBLE * alpha * PI) * cc**alpha
    s     = 1.0_C_DOUBLE / alpha
    IF ( B .LE. 0.0_C_DOUBLE ) RETURN
    x     = B * t0**alpha
    IF ( x .LE. s ) RETURN

    lb = C - x + (s - 1.0_C_DOUBLE) * DLOG(x) + DLOG( x / (x - s + 1.0_C_DOUBLE) ) &
         - DLOG(alpha) - s * DLOG(B)

  END FUNCTION logTailBoundPgt2
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE updateTM(i, tmax, mmax, left_Of_Max, &
                      m, zeroL, zeroR, error, xacc_in, mOld_out, &
                      zeroBoundL_out, zeroBoundR_out)
                  
    ! Update the values of  t  and  m  to the values needed
    ! for the next integration region.
    ! The updated value of  m  that is returned corresponds to the
    ! updated value of  zeroR.
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax
    REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroR
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroL
    INTEGER(C_INT), INTENT(IN)          :: i, mmax
    INTEGER(C_INT), INTENT(INOUT)       :: m
    LOGICAL(C_BOOL), INTENT(INOUT)      :: left_Of_Max, error
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: xacc_in
    INTEGER(C_INT), INTENT(OUT)         :: mOld_out
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroBoundL_out, zeroBoundR_out
      ! The exact bracket [zeroBoundL_out, zeroBoundR_out] used to find
      ! zeroR this call, after improveKZeroBounds refined it. Exposed so
      ! that, if this turns out to be the final pre-acceleration call, the
      ! caller can re-run findExactZeros on this SAME bracket at a tighter
      ! tolerance, rather than guessing a new one.
    
    ! Local vars
    INTEGER(C_INT)        :: mOld
    REAL(KIND=C_DOUBLE)   :: zeroBoundL, zeroBoundR, zeroStartPoint
    REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi

    
    current_y    = Cy(i)
    current_mu   = Cmu(i)
    current_phi  = Cphi(i)

    zeroL = zeroR
    mOld = m
    
    CALL advanceM(m, mmax, mOld, left_Of_Max)

    IF ( left_Of_Max ) THEN
      zeroBoundR = tmax
      zeroBoundL = zeroR
    ELSE 
      zeroBoundL = tmax
      zeroBoundR = zeroBoundL * 20.0_C_DOUBLE
    END IF
    
    zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE

    CALL improveKZeroBounds(m, left_Of_Max, zeroStartPoint, &
                            zeroBoundL, zeroBoundR, error)
    zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE

    CALL findExactZeros(m, zeroBoundL, zeroBoundR, &
                        zeroStartPoint, zeroR, left_Of_Max, error, xacc_in)

    mOld_out = mOld
    zeroBoundL_out = zeroBoundL
    zeroBoundR_out = zeroBoundR
    
  END SUBROUTINE updateTM
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  SUBROUTINE findInitialZeroR(mfirst, left_Of_Max, tmax, &
                              zeroR, error, xacc_in)
  
    USE tweedie_params_mod
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroR
    INTEGER(C_INT), INTENT(IN)          :: mfirst
    LOGICAL(C_BOOL), INTENT(INOUT)      :: left_Of_Max
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax
    LOGICAL(C_BOOL), INTENT(INOUT)      :: error
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: xacc_in

    ! Local vars
    REAL(KIND=C_DOUBLE)                 :: t_Start_Point, zeroBoundL, zeroBoundR
    REAL(KIND=C_DOUBLE)                 :: TMP
    LOGICAL(C_BOOL)                     :: errorHere

    ! Initialisation
    t_Start_Point = 0.0_C_DOUBLE
    zeroBoundL = 0.0_C_DOUBLE
    zeroBoundR = 0.0_C_DOUBLE
    zeroR = 0.0_C_DOUBLE
    TMP = 0.0_C_DOUBLE
    errorHere = .FALSE.
    
    ! Find starting point for the first zero
    IF (left_Of_Max) THEN
      t_Start_Point = PI / current_y  
      zeroBoundL = 0_C_DOUBLE
      zeroBoundR = tmax   ! WAS: t_Start_Point * 2.0_C_DOUBLE
    ELSE
      ! Searching to the right of tmax
      t_Start_Point = tmax + PI / current_y  
      zeroBoundL = tmax
      zeroBoundR = t_Start_Point * 2.0_C_DOUBLE
    END IF

    IF ( (t_Start_Point .GT. zeroBoundR) .OR. (t_Start_Point .LT. zeroBoundL) ) Then
      t_Start_Point = (zeroBoundL + zeroBoundR) / 2.0_c_DOUBLE
    END IF
  
    ! Find the zero
    CALL findExactZeros(mfirst, zeroBoundL, zeroBoundR, t_Start_Point, zeroR, & 
                        left_Of_Max, errorHere, xacc_in)
    ! findExactZeros may change the value of  left_Of_Max

    CALL evaluateImk(zeroR, TMP, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: Imk: integrand zero =", -1, zeroR, 1)
    END IF
    
  END SUBROUTINE findInitialZeroR

END MODULE TweedieIntHelpers

