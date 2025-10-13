SUBROUTINE DFsmallp(funvalue, exitstatus, relerr, verbose) BIND(C, NAME='DFsmallp')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! --- Subroutine Outputs/Inputs ---
  REAL(KIND=C_DOUBLE), INTENT(OUT)     :: funvalue, relerr
  INTEGER(C_INT), INTENT(OUT)   :: exitstatus
  INTEGER(C_INT), INTENT(IN)    :: verbose

  
  ! --- Local Variables ---
  REAL(KIND=8) :: pi, area, zeroL, zeroR, zero, zeroStartPoint
  REAL(KIND=8) :: epsilon, aimrerr, omega
  INTEGER(C_INT) :: m

  ! --- External Function/Subroutine Declarations (CRITICAL: only EXTERNAL needed) ---
  ! DFsmallp only calls DFintegrand, gaussq, and findExactZeros directly.
  EXTERNAL DFintegrand, gaussq, findExactZeros
  REAL(KIND=8) DFintegrand
  
  ! --- Initialization ---
  CpSmall = .TRUE.
  
  pi = 4.0_8 * DATAN(1.0_8)
  exitstatus = 0
  relerr = 0.0_8 
  epsilon = 1.0_8 / 10.0_8**16 ! 1.0d-16
  aimrerr = 1.0_8 / 10.0_8**14 ! 1.0d-14
  
  IF (verbose .EQ. 1) WRITE(*, '(A)') " FOR p <= 2"
  
  ! 1. Define initial integration bounds
  zeroL = 0.0_8
  zeroStartPoint = pi / 2.0_8
  
  ! 2. Find the upper integration limit (m=0 is typical for the first zero)
  m = 0
  ! NOTE: The call to findExactZeros is correct and uses the function pointer.
  CALL findExactZeros(zeroL, 10.0_8, zeroStartPoint, zero)
  zeroR = zero 
  
  IF (verbose .EQ. 1) WRITE(*, '(A, 2F10.8)') "  - Integrate between 0 and ", zeroR
  
  ! 3. Perform the integration
  CALL gaussq( DFintegrand, area, zeroL, zeroR)
  
  ! 4. Calculate final CDF value
  funvalue = (-1.0_8/pi) * area + 0.5_8
  
  ! 5. Set final exit status
  IF (relerr < aimrerr) THEN
    exitstatus = 1
  ELSE
    exitstatus = -10
  END IF
  
  IF (verbose .EQ. 1) THEN
    WRITE(*, '(A, F10.8)') "FINAL CDF VALUE:", funvalue
  END IF
  
END SUBROUTINE DFsmallp
