MODULE R_interfaces
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR
  IMPLICIT NONE

CONTAINS

  SUBROUTINE INTPR(S, N, V, NV)
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: S(*)
    INTEGER, INTENT(IN) :: N
    ! OPTIONAL makes both 2-arg and 4-arg calls valid simultaneously
    INTEGER, INTENT(IN), OPTIONAL :: V
    INTEGER, INTENT(IN), OPTIONAL :: NV
    
    INTEGER(C_INT) :: V_tmp(1), NV_tmp
    
    INTERFACE
       SUBROUTINE intpr_c(S, N, V, NV) BIND(C, NAME="intpr_")
         USE ISO_C_BINDING
         CHARACTER(KIND=C_CHAR) :: S(*)
         INTEGER(C_INT), VALUE :: N
         INTEGER(C_INT), INTENT(IN) :: V(*)
         INTEGER(C_INT), VALUE :: NV
       END SUBROUTINE intpr_c
    END INTERFACE

    IF (PRESENT(V) .AND. PRESENT(NV)) THEN
       V_tmp(1) = INT(V, C_INT)
       NV_tmp = INT(NV, C_INT)
    ELSE
       V_tmp(1) = 0_C_INT
       NV_tmp = 0_C_INT
    END IF

    CALL intpr_c(S, INT(N, C_INT), V_tmp, NV_tmp)
  END SUBROUTINE INTPR

  SUBROUTINE DBLEPR(S, N, V, NV)
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: S(*)
    INTEGER, INTENT(IN) :: N
    REAL(KIND=C_DOUBLE), INTENT(IN), OPTIONAL :: V
    INTEGER, INTENT(IN), OPTIONAL :: NV
    
    REAL(KIND=C_DOUBLE) :: V_tmp(1)
    INTEGER(C_INT) :: NV_tmp
    
    INTERFACE
       SUBROUTINE dblepr_c(S, N, V, NV) BIND(C, NAME="dblepr_")
         USE ISO_C_BINDING
         CHARACTER(KIND=C_CHAR) :: S(*)
         INTEGER(C_INT), VALUE :: N
         REAL(C_DOUBLE), INTENT(IN) :: V(*)
         INTEGER(C_INT), VALUE :: NV
       END SUBROUTINE dblepr_c
    END INTERFACE

    IF (PRESENT(V) .AND. PRESENT(NV)) THEN
       V_tmp(1) = V
       NV_tmp = INT(NV, C_INT)
    ELSE
       V_tmp(1) = 0.0_C_DOUBLE
       NV_tmp = 0_C_INT
    END IF

    CALL dblepr_c(S, INT(N, C_INT), V_tmp, NV_tmp)
  END SUBROUTINE DBLEPR

END MODULE R_interfaces