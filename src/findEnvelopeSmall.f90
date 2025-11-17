SUBROUTINE findEnvelopeSmall(i, t_Start, t_Out)
  ! Finds where the envelope of the integrand is miniscule, in terms of t.
  !
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)           :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)      :: t_Start
  REAL(KIND=C_DOUBLE), INTENT(OUT)     :: t_Out

  REAL(KIND=C_DOUBLE) :: pi
  REAL(KIND=C_DOUBLE) :: tL, tR, aimrerr, t_Start_Pt, t_Small_Envelope
  LOGICAL(C_BOOL)     :: root_Found
  REAL(KIND=C_DOUBLE) :: current_y, current_mu, current_phi

  ! ---------------------------
  ! Abstract interface for the function+derivative
  ! ---------------------------
  INTERFACE
    SUBROUTINE funcd_signature(i_in, t, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)         :: i_in
      REAL(KIND=C_DOUBLE), INTENT(IN)    :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)   :: f, df
    END SUBROUTINE funcd_signature

    ! Interface for rtsafe — expects a Fortran procedure matching funcd_signature
    SUBROUTINE rtsafe(i_in, funcd, x1, x2, xacc, root, root_Found)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i_in
      REAL(C_DOUBLE), INTENT(IN)        :: x1, x2, xacc
      REAL(C_DOUBLE), INTENT(OUT)       :: root
      LOGICAL(C_BOOL), INTENT(OUT)      :: root_Found
      PROCEDURE(funcd_signature)        :: funcd  
    END SUBROUTINE rtsafe

  END INTERFACE

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  ! --- Initialization ---
  aimrerr = 1.0E-09_C_DOUBLE
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

  ! Bounds for root search
  tL = 0.0E-09_C_DOUBLE
  tR = t_Start
  t_Start_Pt = 0.75_C_DOUBLE * t_Start

  ! Call rtsafe: pass the Fortran contained procedure evaluateEnvelope
  CALL rtsafe(i, evaluateEnvelope, tL, tR, aimrerr, t_Small_Envelope, root_Found)

WRITE(*,*) "++++++ SO WE FIND ", t_Small_Envelope, root_Found

  IF ( .NOT.(root_Found)) THEN
    t_Out = t_Start
  ELSE
    t_Out = t_Small_Envelope
  END IF
  
CONTAINS

  SUBROUTINE evaluateEnvelope(i_in, t, f, df)
    ! This must be a subroutine matching funcd_signature (not a function).
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
    IMPLICIT NONE

    INTEGER(C_INT), INTENT(IN)        :: i_in
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    
    REAL(KIND=C_DOUBLE)   :: Rek, Rekd, epsilon

    ! Compute f and df via your existing evaluators:
    CALL evaluateRek(i_in, t, Rek)
    CALL evaluateRekd(i_in, t, Rekd)
    
    epsilon = 1.0E-19_C_DOUBLE
    f = DEXP(Rek)/t - epsilon
    df = ( DEXP( Rek ) / t**2.0_C_DOUBLE) *  &
           (t * Rekd - 1.0_C_DOUBLE)
    
  END SUBROUTINE evaluateEnvelope

END SUBROUTINE findEnvelopeSmall
