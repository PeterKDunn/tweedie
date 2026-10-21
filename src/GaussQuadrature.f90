SUBROUTINE GaussQuadrature(i, a, b, integral_result) 
  ! Perform 512-pt Gaussian-Legendre quadrature to evaluate integrals
  
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE Integrands_MOD, ONLY: Integrands
  USE gaussian_data_mod, ONLY: absc, wts
  
  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: a, b      ! Integration limits
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: integral_result
  INTEGER(C_INT), INTENT(IN)        :: i

  INTEGER(C_INT)                    :: j, npoints
  REAL(KIND=C_DOUBLE)               :: xminus, xplus, fplus, fminus
  
  !!! NOTE: Guassian abscicca and weights in gaussian_data_mod
  
  ! Set up initial parameters
  integral_result = 0.0E0_C_DOUBLE
  npoints = SIZE(absc) ! For 512-pt quadrature (e.g.): symmetry
  
  ! Compute
  DO j = 1, npoints
    ! Adjust abscissae
    xplus  = ( b + a ) / 2.0E0_C_DOUBLE + ( b - a ) / 2.0E0_C_DOUBLE * absc(j) 
    xminus = ( b + a ) / 2.0E0_C_DOUBLE - ( b - a ) / 2.0E0_C_DOUBLE * absc(j) 

    ! Evaluate
    fplus  = Integrands(i, xplus)
    fminus = Integrands(i, xminus)
    integral_result = integral_result + wts(j) * (fplus + fminus)
  END DO
  
  integral_result = integral_result * (b - a) / 2.0E0_C_DOUBLE

END SUBROUTINE GaussQuadrature


SUBROUTINE GaussQuadratureAbs(i, a, b, integral_result, abs_result)
  ! As GaussQuadrature, but also returns the same rule applied to |integrand|,
  ! which sets the size of the rounding error in integral_result

  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE Integrands_MOD, ONLY: Integrands
  USE gaussian_data_mod, ONLY: absc, wts

  IMPLICIT NONE

  REAL(KIND=C_DOUBLE), INTENT(IN)   :: a, b
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: integral_result, abs_result
  INTEGER(C_INT), INTENT(IN)        :: i

  INTEGER(C_INT)                    :: j
  REAL(KIND=C_DOUBLE)               :: xminus, xplus, fplus, fminus

  integral_result = 0.0E0_C_DOUBLE
  abs_result      = 0.0E0_C_DOUBLE
  DO j = 1, SIZE(absc)
    xplus  = ( b + a ) / 2.0E0_C_DOUBLE + ( b - a ) / 2.0E0_C_DOUBLE * absc(j)
    xminus = ( b + a ) / 2.0E0_C_DOUBLE - ( b - a ) / 2.0E0_C_DOUBLE * absc(j)
    fplus  = Integrands(i, xplus)
    fminus = Integrands(i, xminus)
    integral_result = integral_result + wts(j) * (fplus + fminus)
    abs_result      = abs_result      + wts(j) * (DABS(fplus) + DABS(fminus))
  END DO
  integral_result = integral_result * (b - a) / 2.0E0_C_DOUBLE
  abs_result      = abs_result      * (b - a) / 2.0E0_C_DOUBLE

END SUBROUTINE GaussQuadratureAbs


SUBROUTINE GaussAdaptive(i, a, b, absScale, aimrerr, area, ok)
  ! Integrate over [a, b] with the 512-point Gauss rule, subdividing where
  ! needed (globally adaptive, as in QUADPACK's QAG).
  !
  ! WHY: an integration region runs between consecutive zeros of the
  ! integrand, found from the levels  Im k(t) = m pi. The zero-finding assumes
  ! that Im k(t) passes each level once on each side of its maximum. For p near
  ! 1 that is false: alpha = (2-p)/(1-p) is large and negative, so
  ! (1 + i x)^alpha winds round the origin many times while it decays, and
  ! Im k(t) oscillates about zero before the  -t y  term takes over. A single
  ! region can then hold many oscillations of the integrand, and one Gauss
  ! rule over it is wrong. (E.g. ptweedie_inversion(0.001, mu = 5, phi = 10,
  ! power = 1.01): the region [7.9, 15708] (scaled to mu = 1) has all its
  ! oscillations in [7.9, 200], where the rule has 36 of its 512 nodes; F was
  ! in error by 2e-5, with no warning.)
  !
  ! HOW: each interval's value is the sum of the rule over its two halves; its
  ! error estimate is the difference from the rule over the whole interval.
  ! The interval with the largest estimated error is halved, until the total
  ! estimated error is below  max(aimrerr, 1e-13) * scale, where scale is the larger of
  ! absScale (sum of |areas| integrated before this region) and the integral
  ! of |integrand| over [a, b]. The second sets the rounding error, so the
  ! test does not chase rounding noise when the area cancels to near zero.
  ! For a well-behaved region this costs 3 rules (not 1) and changes nothing
  ! beyond rounding. Used only for 1 < p < 4/3 (see below).
  !
  ! ok is .FALSE. if, after MAXINT intervals, the estimated error still
  ! exceeds 1e-8 * scale (a smaller remainder is rounding noise).

  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE tweedie_params_mod, ONLY: CpSmall, Cp

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: a, b, absScale, aimrerr
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: area
  LOGICAL(C_BOOL), INTENT(OUT)      :: ok

  INTEGER, PARAMETER                :: MAXINT = 40
  REAL(KIND=C_DOUBLE), PARAMETER    :: QTOL_FLOOR = 1.0E-13_C_DOUBLE
  REAL(KIND=C_DOUBLE), PARAMETER    :: FLAG_REL   = 1.0E-08_C_DOUBLE
    ! The two estimates compared (one rule; two half rules) each carry
    ! rounding error of order 1e-16 * sqrt(512) relative to the integral of
    ! |integrand|, so a target below about 1e-13 would chase rounding noise
  REAL(KIND=C_DOUBLE)               :: lo(MAXINT), hi(MAXINT), val(MAXINT), err(MAXINT)
  REAL(KIND=C_DOUBLE)               :: halfL(MAXINT), halfR(MAXINT), absv(MAXINT)
  REAL(KIND=C_DOUBLE)               :: whole, mid, tol, scale
  INTEGER                           :: n, j, k

  INTERFACE
    SUBROUTINE GaussQuadrature(i, a, b, area)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(KIND=C_DOUBLE), INTENT(OUT)      :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)       :: a, b
      INTEGER(C_INT), INTENT(IN)            :: i
    END SUBROUTINE GaussQuadrature
    SUBROUTINE GaussQuadratureAbs(i, a, b, area, absArea)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(KIND=C_DOUBLE), INTENT(OUT)      :: area, absArea
      REAL(KIND=C_DOUBLE), INTENT(IN)       :: a, b
      INTEGER(C_INT), INTENT(IN)            :: i
    END SUBROUTINE GaussQuadratureAbs
  END INTERFACE

  ok = .TRUE.
  IF (b .LE. a) THEN
    area = 0.0_C_DOUBLE
    RETURN
  END IF

  ! Only p < 4/3 needs this. Write  Im k(t) = lambda g(x) - t y,  x = |c| t,
  ! theta = atan(x), A = |alpha| = (2-p)/(p-1); for 1 < p < 2,
  !     g = cos(theta)^A sin(A theta),
  !     dg/dx = A cos(theta)^(A+1) cos((A+1) theta).
  ! If A <= 2 (p >= 4/3), (A+1) theta < 3 pi/2, so dg/dx is decreasing while
  ! positive and stays negative after the peak: Im k(t) then has a single
  ! maximum, each level m pi is crossed once on each side, and the regions
  ! between the zeros found are well behaved. For p > 2, alpha is in (0, 1)
  ! and the same holds. In those cases: one rule, as before.
  IF ( (.NOT. CpSmall) .OR. (Cp .GE. 4.0_C_DOUBLE / 3.0_C_DOUBLE) ) THEN
    CALL GaussQuadrature(i, a, b, area)
    RETURN
  END IF

  n = 1
  lo(1) = a
  hi(1) = b
  CALL GaussQuadrature(i, a, b, whole)
  CALL splitOne(1, whole)

  DO
    area  = SUM(val(1:n))
    scale = MAX(absScale, SUM(absv(1:n)))
    tol   = MAX(aimrerr, QTOL_FLOOR) * scale
    IF (SUM(err(1:n)) .LE. tol) EXIT
    IF (n .GE. MAXINT) THEN
      ! Out of intervals. What error remains is usually rounding in the
      ! integrand itself (e.g. sin of a large argument), which no subdivision
      ! removes; flag only an error too large for that
      ok = ( SUM(err(1:n)) .LE. FLAG_REL * scale )
      EXIT
    END IF
    ! Halve the interval with the largest estimated error
    j     = MAXLOC(err(1:n), DIM = 1)
    n     = n + 1
    k     = n
    mid   = 0.5_C_DOUBLE * (lo(j) + hi(j))
    lo(k) = mid
    hi(k) = hi(j)
    hi(j) = mid
    whole = halfR(j)
    CALL splitOne(k, whole)       ! right half (new interval k)
    whole = halfL(j)
    CALL splitOne(j, whole)       ! left half (replaces interval j)
  END DO

CONTAINS

  SUBROUTINE splitOne(j, wholeJ)
    ! Value and error estimate of interval j, given the rule over all of it
    INTEGER, INTENT(IN)             :: j
    REAL(KIND=C_DOUBLE), INTENT(IN) :: wholeJ
    REAL(KIND=C_DOUBLE)             :: m, aL, aR
    m = 0.5_C_DOUBLE * (lo(j) + hi(j))
    CALL GaussQuadratureAbs(i, lo(j), m, halfL(j), aL)
    CALL GaussQuadratureAbs(i, m, hi(j), halfR(j), aR)
    val(j)  = halfL(j) + halfR(j)
    err(j)  = DABS(val(j) - wholeJ)
    absv(j) = aL + aR             ! integral of |integrand| over interval j
  END SUBROUTINE splitOne

END SUBROUTINE GaussAdaptive
