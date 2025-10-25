MODULE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE , C_BOOL
  
  IMPLICIT NONE
  SAVE
  
  ! Global Parameters
  REAL(KIND=C_DOUBLE), ALLOCATABLE :: Cmu(:), Cphi(:), Cy(:) 
  REAL(KIND=C_DOUBLE)   :: Cp
  LOGICAL(C_BOOL)       :: CpSmall, Cverbose
  INTEGER(C_INT)        :: CN, current_i
  
END MODULE tweedie_params_mod

