SUBROUTINE findExactZeros(i, m, xL, xR, xStart, xZero) 
  USE tweedie_params_mod, ONLY: current_i
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  IMPLICIT NONE

  ! --- Dummy Arguments ---
  INTEGER(C_INT), INTENT(IN)        :: i, m             ! i: which value; m: the m*pi value to be solving for
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: xL, xR, xStart   ! L and R bounds, and starting point
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: xZero
  
  ! --- Local Variables
  REAL(KIND=C_DOUBLE) :: xacc ! Accuracy value


  INTERFACE
    ! The signature required by rtnewton (3 arguments: x, f, df)
    SUBROUTINE funcd_signature(i, x, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE 
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: x
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    END SUBROUTINE funcd_signature

    SUBROUTINE rtnewton(i, funcd, x1, x2, xstart, xacc, root)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE

      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: x1, x2, xstart, xacc
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root
      
      PROCEDURE(funcd_signature) :: funcd
      
    END SUBROUTINE rtnewton

    ! 2. The Function being Solved (findImkM)
    ! CRITICAL FIX: Signature reduced to 3 arguments (t, f, df) to match funcd_signature
    SUBROUTINE findImkM(i, t, f, df, m)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      
      INTEGER(C_INT)                    :: i, m
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    END SUBROUTINE findImkM

  END INTERFACE
  ! --- END INTERFACES ---
  
  
!  WRITE(*,*) "findExactZeros: INSIDE"
  ! Set the accuracy
  xacc = 1.0D-12
  
  ! CRITICAL STEP: Set the global context variable 'current_i' before calling the root finder.
  current_i = i
!  WRITE(*,*) "findExactZeros: ABOUT to call RTNEWTON"
  
  CALL rtnewton(i, findImkM_wrapper, xL, xR, xStart, xacc, xZero)
!  WRITE(*,*) "findExactZeros: Called RTNEWTON, about to exit findExactZeros"
  
  CONTAINS 
  
    ! Define the wrapper subroutine
    SUBROUTINE findImkM_wrapper(i, x, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      ! Has access to bariable  m  from the containing function/outer scope
      
      INTEGER(C_INT), INTENT(IN)  :: i
      REAL(C_DOUBLE), INTENT(IN)  :: x
      REAL(C_DOUBLE), INTENT(OUT) :: f, df
     
      CALL findImkM(i, x, f, df, m) ! Pass the captured 'm' value
    END SUBROUTINE findImkM_wrapper

  
END SUBROUTINE findExactZeros
