SUBROUTINE findExactZeros(i, m, xL, xR, xStart, xZero) 

  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  IMPLICIT NONE

  ! --- Dummy Arguments ---
  INTEGER(C_INT), INTENT(IN)        :: i, m             ! i: which value; m: the m*pi value to be solving for
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: xL, xR, xStart   ! L and R bounds, and starting point
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: xZero
  
  ! --- Local Variables
  REAL(KIND=C_DOUBLE) :: xacc ! Accuracy value


  INTERFACE
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


    SUBROUTINE findImkM(i, t, f, df, m)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i, m
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    END SUBROUTINE findImkM

  END INTERFACE

!write(*,*) "** IN fundxacteros!!"
  
  ! Set the accuracy
  xacc = 1.0E-13_C_DOUBLE
!write(*,*) "** IN findexactzeros1, ",  xL, xR

  
  CALL rtnewton(i, findImkM_wrapper, xL, xR, xStart, xacc, xZero)
!write(*,*) "** IN findexactzeros2, ", xZero

  CONTAINS 
  
    ! Define the wrapper subroutine
    SUBROUTINE findImkM_wrapper(i, x, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      ! Has access to variable  m  from the containing function/outer scope
      
      INTEGER(C_INT), INTENT(IN)  :: i
      REAL(C_DOUBLE), INTENT(IN)  :: x
      REAL(C_DOUBLE), INTENT(OUT) :: f, df
     
      CALL findImkM(i, x, f, df, m) ! Pass the captured 'm' value]
!      write(*,*) "** IN findImk_wrapper, x, f and df:", x, f, df
      
    END SUBROUTINE findImkM_wrapper

  
END SUBROUTINE findExactZeros
