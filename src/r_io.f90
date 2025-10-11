MODULE R_IO_MOD
  USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE Rprintf(fmt, str) BIND(C, NAME='Rprintf')
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: fmt
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: str
    END SUBROUTINE Rprintf
  END INTERFACE
  
  ! You can optionally define a helper subroutine here for easier calls
  CONTAINS
    SUBROUTINE R_PRINT_MSG(message)
      CHARACTER(*), INTENT(IN) :: message
      CHARACTER(KIND=C_CHAR, LEN=LEN(message)+1) :: c_msg
      
      c_msg = TRIM(message) // C_NULL_CHAR
      CALL Rprintf(c_msg, C_NULL_CHAR)
    END SUBROUTINE R_PRINT_MSG

END MODULE R_IO_MOD
