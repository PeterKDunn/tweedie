MODULE R_interfaces
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR
  IMPLICIT NONE

CONTAINS

  ! One routine to rule them all - No generic interface to cause ambiguity
  SUBROUTINE DBLEPR(S, N, V, NV)
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: S(*)
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: V
    INTEGER, INTENT(IN) :: NV
    
    ! Internal C-bound interface
    INTERFACE
       SUBROUTINE dblepr_c(S, N, V, NV) BIND(C, NAME="dblepr_")
         USE ISO_C_BINDING
         CHARACTER(KIND=C_CHAR) :: S(*)
         INTEGER(C_INT), VALUE :: N
         REAL(C_DOUBLE), INTENT(IN) :: V(*)
         INTEGER(C_INT), VALUE :: NV
       END SUBROUTINE dblepr_c
    END INTERFACE

    ! Force conversion to the types R's C code expects
    CALL dblepr_c(S, INT(N, C_INT), [REAL(V, C_DOUBLE)], INT(NV, C_INT))
  END SUBROUTINE DBLEPR

  SUBROUTINE INTPR(S, N, V, NV)
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: S(*)
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: V
    INTEGER, INTENT(IN) :: NV
    
    INTERFACE
       SUBROUTINE intpr_c(S, N, V, NV) BIND(C, NAME="intpr_")
         USE ISO_C_BINDING
         CHARACTER(KIND=C_CHAR) :: S(*)
         INTEGER(C_INT), VALUE :: N
         INTEGER(C_INT), INTENT(IN) :: V(*)
         INTEGER(C_INT), VALUE :: NV
       END SUBROUTINE intpr_c
    END INTERFACE

    CALL intpr_c(S, INT(N, C_INT), [INT(V, C_INT)], INT(NV, C_INT))
  END SUBROUTINE INTPR

END MODULE R_interfaces