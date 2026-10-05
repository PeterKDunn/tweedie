
MODULE tweedie_params_mod
  ! Save commonly-used parameters, for access throughout functions and subroutines
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE , C_BOOL
  
  IMPLICIT NONE
  SAVE
  
  ! Global Parameters
  REAL(KIND=C_DOUBLE), ALLOCATABLE  :: Cmu(:), Cphi(:), Cy(:)     ! mu, phi, y (the C means COMMON) as vectors
  REAL(KIND=C_DOUBLE)               :: Cp                         ! p
  LOGICAL(C_BOOL)                   :: CpSmall                    ! Flag to indicate if p is small (1 < p < 2)
  LOGICAL(C_BOOL)                   :: Cverbose                   ! Flag to indicate if verbose output is requested
  LOGICAL(C_BOOL)                   :: Cpdf                       ! Flag to indicate if the PDF (.TRUE.) or CDF (.FALSE.) is requested
  LOGICAL(C_BOOL)                   :: Ctail                      ! Flag to indciate if teh lower or upper tail is needed for the PDF:
                                                                  !    .FALSE. = lower tail F(y);  .TRUE. = upper tail S(y) = 1 = F(y)
                                                                  !    Has no effect if the PDF is requested
  INTEGER(C_INT)                    :: CN                         ! The number of points for which the evaluation is needed
  

  REAL(KIND=C_DOUBLE) :: current_y, current_mu, current_phi        
  INTEGER(C_INT)      :: m_shared

  REAL(KIND=C_DOUBLE), PARAMETER :: PI = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

CONTAINS

  PURE SUBROUTINE zPow(x, e, re, im)
    ! Real and imaginary parts of  z^e,  where  z = 1 + i x  (for x real).
    !
    ! Equivalently, with  omega = atan(x):
    !   re = cos(e*omega) / cos(omega)^e,   
    !   im = sin(e*omega) / cos(omega)^e,
    ! which is how Re k, Im k and their derivatives were written before.
    ! Computing them that way loses all accuracy for large |x|:
    !  - once |x| > ~1e16, DATAN(x) rounds to exactly +/- pi/2, so DCOS(omega)
    !    returns ~6e-17 instead of ~1/|x|;
    !  - cos(e*omega) or sin(e*omega) can be a tiny value obtained by
    !    cancellation (e.g. e = -1, where cos(e*omega) ~ 1/|x|).
    ! Here, for |x| > 1, omega = s(pi/2 - delta) with delta = atan(1/|x|)
    ! (accurate), and the trig is expanded with angle-sum formulas using
    ! exact sin(pi r) / cos(pi r), so e.g. cos(pi/2) is exactly 0.
    ! |z|^e is computed on the log scale.

    REAL(KIND=C_DOUBLE), INTENT(IN)   :: x, e
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: re, im
    REAL(KIND=C_DOUBLE)               :: logMod, w, delta, s, A, B, modE, cA, sA

    IF (x .EQ. 0.0_C_DOUBLE) THEN
      re = 1.0_C_DOUBLE
      im = 0.0_C_DOUBLE
      RETURN
    END IF

    IF (DABS(x) .LE. 1.0_C_DOUBLE) THEN
      logMod = 0.5_C_DOUBLE * DLOG(1.0_C_DOUBLE + x * x)
      modE   = DEXP(e * logMod)
      w      = DATAN(x)
      re     = modE * DCOS(e * w)
      im     = modE * DSIN(e * w)
    ELSE
      logMod = DLOG(DABS(x)) + 0.5_C_DOUBLE * DLOG(1.0_C_DOUBLE + 1.0_C_DOUBLE / (x * x))
      modE   = DEXP(e * logMod)
      s      = SIGN(1.0_C_DOUBLE, x)
      delta  = DATAN(1.0_C_DOUBLE / DABS(x))
      ! e*omega = pi*A - B,  with  A = s*e/2  and  B = s*e*delta
      A  = s * e / 2.0_C_DOUBLE
      B  = s * e * delta
      cA = cosPi(A)
      sA = sinPi(A)
      re = modE * ( cA * DCOS(B) + sA * DSIN(B) )
      im = modE * ( sA * DCOS(B) - cA * DSIN(B) )
    END IF

  END SUBROUTINE zPow


  PURE FUNCTION sinPi(r) RESULT(v)
    ! sin(pi * r), exact at multiples of 1/2
    REAL(KIND=C_DOUBLE), INTENT(IN) :: r
    REAL(KIND=C_DOUBLE)             :: v, a, sgn
    a   = r - 2.0_C_DOUBLE * ANINT(r / 2.0_C_DOUBLE)     ! now in [-1, 1]
    sgn = SIGN(1.0_C_DOUBLE, a)
    a   = DABS(a)
    IF (a .LE. 0.25_C_DOUBLE) THEN
      v = DSIN(PI * a)
    ELSE IF (a .LE. 0.75_C_DOUBLE) THEN
      v = DCOS(PI * (a - 0.5_C_DOUBLE))
    ELSE
      v = DSIN(PI * (1.0_C_DOUBLE - a))
    END IF
    v = sgn * v
  END FUNCTION sinPi


  PURE FUNCTION cosPi(r) RESULT(v)
    ! cos(pi * r), exact at multiples of 1/2
    REAL(KIND=C_DOUBLE), INTENT(IN) :: r
    REAL(KIND=C_DOUBLE)             :: v, a
    a = DABS(r - 2.0_C_DOUBLE * ANINT(r / 2.0_C_DOUBLE))  ! now in [0, 1]
    IF (a .LE. 0.25_C_DOUBLE) THEN
      v = DCOS(PI * a)
    ELSE IF (a .LE. 0.75_C_DOUBLE) THEN
      v = -DSIN(PI * (a - 0.5_C_DOUBLE))
    ELSE
      v = -DCOS(PI * (1.0_C_DOUBLE - a))
    END IF
  END FUNCTION cosPi

END MODULE tweedie_params_mod

