!======================================================================
! File: interface.f90
! Purpose: Expose the main twcdf subroutine to R using BIND(C)
!======================================================================
SUBROUTINE twcdf(N, p, phi, y, mu, funvalue, exitstatus, relerr, its) BIND(C, name='twcdf')
  
  ! This file must use the module where your data structures are defined.
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  ! --- Arguments for the R/C Interface ---
  ! Note: R passes all arguments by reference, which naturally maps to Fortran.
  
  INTEGER(C_INT), INTENT(IN)      :: N          
  REAL(KIND=C_DOUBLE), INTENT(IN)      :: p          
  REAL(KIND=C_DOUBLE), INTENT(IN)      :: phi(*)     ! Assumed-size array for R interface
  REAL(KIND=C_DOUBLE), INTENT(IN)      :: y(*)       
  REAL(KIND=C_DOUBLE), INTENT(IN)      :: mu(*)      
  REAL(KIND=C_DOUBLE), INTENT(OUT)     :: funvalue(*) 
  
  INTEGER(C_INT), INTENT(OUT)   :: exitstatus 
  REAL(KIND=C_DOUBLE), INTENT(OUT)   :: relerr     
  INTEGER(C_INT), INTENT(OUT)   :: its        
  
  ! ------------------------------------------------------------------
  
  ! Call the original Fortran subroutine.
  ! We use the original 'twcdf' name here. The compiler will handle 
  ! the difference between the BIND(C) entry point and the internal call.
  ! The argument lists must match exactly.
  CALL twcdf_internal(N, p, phi, y, mu, funvalue, exitstatus, relerr, its)
  
END SUBROUTINE twcdf


!======================================================================
! RENAME: Your original twcdf subroutine must be renamed to twcdf_internal
!         (or similar) and moved to a separate file (e.g., twcdf_main.f90).
!         This prevents the compiler from seeing two routines named 'twcdf'.
!======================================================================
! You would NOT put the internal code here, but it's shown for context:
!
! SUBROUTINE twcdf_internal(N, p, phi, y, mu, funvalue, exitstatus, relerr, its)
!   USE tweedie_params_mod
!   IMPLICIT NONE
!   ... (All your local declarations and the main loop logic go here) ...
!   ...
! END SUBROUTINE twcdf_internal