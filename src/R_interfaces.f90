MODULE R_interfaces
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR
  IMPLICIT NONE

  INTERFACE DBLEPR
     MODULE PROCEDURE DBLEPR_SCALAR
     MODULE PROCEDURE DBLEPR_ARRAY
  END INTERFACE DBLEPR

CONTAINS

  SUBROUTINE DBLEPR_SCALAR(S, N, V, NV)
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: S
    INTEGER(C_INT), VALUE :: N
    REAL(C_DOUBLE), INTENT(IN) :: V
    INTEGER(C_INT), VALUE :: NV
    
    ! Internal interface to R's actual printing routine
    INTERFACE
       SUBROUTINE dblepr(S, N, V, NV) BIND(C, NAME="dblepr_") 
         USE ISO_C_BINDING
         CHARACTER(KIND=C_CHAR) :: S(*)
         INTEGER(C_INT), VALUE :: N
         REAL(C_DOUBLE), INTENT(IN) :: V
         INTEGER(C_INT), VALUE :: NV
       END SUBROUTINE dblepr
    END INTERFACE
    
    CALL dblepr(S, N, V, NV)
  END SUBROUTINE DBLEPR_SCALAR

  SUBROUTINE DBLEPR_ARRAY(S, N, V, NV)
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: S
    INTEGER(C_INT), VALUE :: N
    REAL(C_DOUBLE), INTENT(IN) :: V(*)
    INTEGER(C_INT), VALUE :: NV
    
    INTERFACE
       SUBROUTINE dblepr(S, N, V, NV) BIND(C, NAME="dblepr_") 
         USE ISO_C_BINDING
         CHARACTER(KIND=C_CHAR) :: S(*)
         INTEGER(C_INT), VALUE :: N
         REAL(C_DOUBLE), INTENT(IN) :: V(*)
         INTEGER(C_INT), VALUE :: NV
       END SUBROUTINE dblepr
    END INTERFACE
    
    CALL dblepr(S, N, V, NV)
  END SUBROUTINE DBLEPR_ARRAY
  
  SUBROUTINE INTPR(S, N)
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: S(*)
    INTEGER(C_INT), VALUE :: N
    
    INTERFACE
       SUBROUTINE intpr_c(S, N, V, NV) BIND(C, NAME="intpr_")
         USE ISO_C_BINDING
         CHARACTER(KIND=C_CHAR) :: S(*)
         INTEGER(C_INT), VALUE :: N
         INTEGER(C_INT) :: V  ! Dummy for intpr
         INTEGER(C_INT), VALUE :: NV
       END SUBROUTINE intpr_c
    END INTERFACE
  
    ! R's intpr(label, n_label, vector, n_vector)
    ! We pass 0 and 0 for the vector parts when just printing a string
    CALL intpr_c(S, N, 0_C_INT, 0_C_INT)
  END SUBROUTINE INTPR

END MODULE R_interfaces