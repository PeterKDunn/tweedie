
SUBROUTINE twcdf_internal(N, p, phi, y, mu, funvalue, exitstatus, relerr, its) BIND(C, NAME='twcdf_internal')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! --- Subroutine Outputs ---
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus, its

  INTEGER(C_INT), INTENT(IN)        :: N
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: p
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: phi(N)
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: y(N)
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: mu(N)
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalue(N)
  
  
  ! --- Local Variables ---
  INTEGER             :: i
  REAL(KIND=8)        :: lambda, aimrerr
  INTEGER             :: verbose
  
  ! --- External Subroutines (Must be updated to accept all necessary parameters) ---
  EXTERNAL findLambda, DFsmallp, DFbigp

  ! --- Initialization: These are then common in 00tweedie_params.F90
  Cp = p
  Cy = y
  Cmu = mu
  Cphi = phi
  CN = N

  verbose = 1
  exitstatus = 1
  relerr = 0.0_8
  its = 0
  aimrerr = 1.0_8 / 10.0_8**10 ! 1.0d-10

  
  ! --- Determine case: psmall = TRUE means 1 < p < 2 ---
  CpSmall = .FALSE.
  IF ( (p > 1.0_8 ) .AND. (p < 2.0_8) ) CpSmall = .TRUE.
  

  ! --- Compute lambda (pSmall case only) ---
  ! NOTE: findLambda must be updated to accept p, mu, phi, and potentially y
  CALL findLambda(lambda, p, mu, phi)
  
  ! --- SPECIAL CASE: if y < 0, return 0 ---
!  IF ( y < 0.0_8 ) THEN
!    funvalue = 0.0_8
!    RETURN
!  END IF
  
  ! --- SPECIAL CASE: if 1 < p < 2, Pr(Y = 0) = exp( -lambda ) ---
!  IF ( CpSmall .AND. (y == 0.0_8 ) ) THEN
!    funvalue = EXP( -lambda )
!    RETURN
!  END IF
  
  
  ! If you desire runtime logging, uncomment the following block
  ! WRITE(*, '(A, F10.8)') "** Computing for y: ", y
  ! WRITE(*, '(A, F10.8)') "** mu: ", mu
  ! WRITE(*, '(A, F10.8)') "** p: ", p
  ! WRITE(*, '(A, F10.8)') "** phi: ", phi
  
  ! --- Call appropriate distribution function routine ---
  DO i = 1, N
    IF ( CpSmall ) THEN
      CALL DFsmallp(i, funvalue, exitstatus, relerr, verbose)
    ELSE
      CALL DFbigp(i, funvalue, exitstatus, relerr, verbose)
    END IF
  END DO
  ! --- Final machine accuracy fixes (currently commented out as in original) ---
  ! IF (psmall) THEN
  !   IF (funvalue < EXP(-lambda) ) funvalue = EXP(-lambda)
  ! ELSE
  !   IF (funvalue < 0.0_8) funvalue = 0.0_8
  ! ENDIF

  ! WRITE(*, '(A)') "IN twcdf: funvalue, exitstatus, relerr"
  ! WRITE(*, '(F10.8, I0, F10.8)') funvalue, exitstatus, relerr

END SUBROUTINE twcdf_internal
