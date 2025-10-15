SUBROUTINE findExactZeros(i, xL, xR, xStart, xZero) 
  USE tweedie_params_mod, ONLY: current_i
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  ! --- Dummy Arguments ---
  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: xL, xR, xStart
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: xZero
  
  ! --- Local Variables
  REAL(KIND=C_DOUBLE) :: xacc ! Accuracy value


  INTERFACE
    ! The signature required by rtnewton (3 arguments: x, f, df)
    SUBROUTINE funcd_signature(x, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE 
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: x
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
    END SUBROUTINE funcd_signature

    FUNCTION rtnewton(funcd, x1, x2, xstart, xacc)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE) :: rtnewton
      
      PROCEDURE(funcd_signature) :: funcd
      
      REAL(KIND=C_DOUBLE), INTENT(IN) :: x1, x2, xstart, xacc
    END FUNCTION rtnewton

    ! 2. The Function being Solved (findImkM)
    ! CRITICAL FIX: Signature reduced to 3 arguments (t, f, df) to match funcd_signature
    SUBROUTINE findImkM(t, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      
      REAL(KIND=C_DOUBLE), INTENT(IN)  :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: f, df
    END SUBROUTINE findImkM

  END INTERFACE
  ! --- END INTERFACES ---
  
  
  ! Set the accuracy
  xacc = 1.0D-12
  
  ! CRITICAL STEP: Set the global context variable 'current_i' before calling the root finder.
  current_i = i
  
  xZero = rtnewton(findImkM, xL, xR, xStart, xacc)
  
END SUBROUTINE findExactZeros
