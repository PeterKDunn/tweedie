MODULE tweedie_params_mod

  IMPLICIT NONE
  SAVE
  
  ! Global Parameters
  DOUBLE PRECISION :: Cp
  DOUBLE PRECISION, ALLOCATABLE :: Cmu(:), Cphi(:), Cy(:) 
  LOGICAL :: CpSmall
  INTEGER :: CN, current_i
  
END MODULE tweedie_params_mod

