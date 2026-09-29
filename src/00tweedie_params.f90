
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

END MODULE tweedie_params_mod

