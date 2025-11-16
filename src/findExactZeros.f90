SUBROUTINE findExactZeros(i, m, mmax, tmax, tL, tR, tStart, tZero, leftOfMax) 
  ! Find the exact zeros of the integrand
  
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE Calcs_Solvers
  USE Calcs_Imag
  
  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)        :: i, m, mmax
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: tL, tR, tStart, tmax
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: tZero
  LOGICAL(C_BOOL), INTENT(IN)       :: leftOfMax
  
  REAL(KIND=C_DOUBLE)   :: xacc, fL, fR, dfL, dfR, tstart_update, tMid
  REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  ! Set the accuracy
  xacc = 1.0E-13_C_DOUBLE

  ! Find mmax. which depends on whether we are working with the PDF or the CDF.
  ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
  ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m * pi/y;
  !       the CDF has integrand zeros at Im k(t) =        m * pi/y.

  ! Ensure the bounds actually bound the zero
  CALL evaluateImkM(i, tL, fL, dfL, m)
  CALL evaluateImkM(i, tR, fR, dfR, m)

  IF ( (fL * fR) .GT. 0.0_C_DOUBLE ) THEN
    ! Then bounds do not bound the zero
     tMid = (tL + tR) / 2.0_C_DOUBLE
 
    CALL improveKZeroBounds(i, m, leftOfMax, mmax, tmax, tMid, tL, tR)
    IF (Cverbose) CALL DBLEPR("Bounds do not bracket the zero (findExactZeros)", -1, fR, 1)
  END IF

  ! For robustness, use rtsafe when the  distance between zeros 
  ! is expected to be small (i.e., in the tail).
  IF ( m .LE. -3 ) THEN ! Use rtsafe whenever m is large and negative
    CALL rtsafe(i, evaluateImkM_wrapper, tL, tR, xacc, tZero)
  ELSE IF ( (Cpsmall) .AND. (current_y .LT. current_mu) ) THEN
    ! When small p and small y, fight harder for good starting bounds
    CALL improveKZeroBounds(i, m, leftOfMax, mmax, tmax, tStart, tL, tR)
    CALL rtsafe(i, evaluateImkM_wrapper, tL, tR, xacc, tZero)
  ELSE
    ! Default to rtnewton for "easy" cases (e.g., initial zeros)
    tstart_update = (tL + tR) / 2.0_C_DOUBLE
    CALL rtnewton(i, evaluateImkM_wrapper, tstart_update, xacc, tZero)
  END IF

  CONTAINS 
  
    ! Define the wrapper subroutine
    SUBROUTINE evaluateImkM_wrapper(i, x, f, df)
      ! Evaluate  Im k(t) - m * pi (CDF) - pi/2  or  Im k(t) - m * pi (CDF)

      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      ! Has access to variable  m  from the containing function/outer scope
      
      INTEGER(C_INT), INTENT(IN)  :: i
      REAL(C_DOUBLE), INTENT(IN)  :: x
      REAL(C_DOUBLE), INTENT(OUT) :: f, df
     
      CALL evaluateImkM(i, x, f, df, m) ! Pass the captured 'm' value]

    END SUBROUTINE evaluateImkM_wrapper

  
END SUBROUTINE findExactZeros
