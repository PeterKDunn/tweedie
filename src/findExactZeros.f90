SUBROUTINE findExactZeros(i, m, tL, tR, tStart, tZero, leftOfMax) 
  USE tweedie_params_mod

  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  IMPLICIT NONE

  ! --- Dummy Arguments ---
  INTEGER(C_INT), INTENT(IN)        :: i, m, leftOfMax  ! i: which value; m: the m*pi value to be solving for
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: tL, tR, tStart   ! L and R bounds, and starting point
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: tZero
  
  ! --- Local Variables
  REAL(KIND=C_DOUBLE) :: xacc, fL, fR, dfL, dfR, current_y, current_mu, current_phi


  INTERFACE
    SUBROUTINE funcd_signature(i, x, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE 
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: x
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    END SUBROUTINE funcd_signature


    SUBROUTINE rtnewton(i, funcd, xstart,xacc, root)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: xstart, xacc
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root
      
      PROCEDURE(funcd_signature) :: funcd
      
    END SUBROUTINE rtnewton


    SUBROUTINE rtsafe(i, funcd, xstart, x1, x2, xacc, root)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: x1, x2, xstart, xacc
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root
      
      PROCEDURE(funcd_signature) :: funcd
      
    END SUBROUTINE rtsafe


    SUBROUTINE findImkM(i, t, f, df, m)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i, m
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    END SUBROUTINE findImkM

  END INTERFACE


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  ! Set the accuracy
  xacc = 1.0E-13_C_DOUBLE

  ! Find mmax. which depends on whether we are working with the PDF or the CDF.
  ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
  ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m pi/y;
  !       the CDF has integrand zeros at Im k(t) =        m pi/y.

  ! Ensure the bounds actually bound the zero
  CALL findImkM(i, tL, fL, dfL, m)
  CALL findImkM(i, tR, fR, dfR, m)
  
  IF ( (fl * fR) .GT. 0.0_C_DOUBLE ) THEN
    ! Then bounds do not actually bound teh ero, so try harder
  END IF

  IF ( (Cpsmall) .AND. (current_y .LT. current_mu) ) THEN
    ! Trickier for small p, small y
    CALL improveKZeroBounds(i, m, leftOfMax, tStart, tL, tR)
    CALL rtsafe(i, findImkM_wrapper, tStart, tL, tR, xacc, tZero)
  ELSE
    ! Newton's method should work fine for big p
    CALL rtnewton(i, findImkM_wrapper, tstart, xacc, tZero)
  END IF

  CONTAINS 
  
    ! Define the wrapper subroutine
    SUBROUTINE findImkM_wrapper(i, x, f, df)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      ! Has access to variable  m  from the containing function/outer scope
      
      INTEGER(C_INT), INTENT(IN)  :: i
      REAL(C_DOUBLE), INTENT(IN)  :: x
      REAL(C_DOUBLE), INTENT(OUT) :: f, df
     
      CALL findImkM(i, x, f, df, m) ! Pass the captured 'm' value]

    END SUBROUTINE findImkM_wrapper

  
END SUBROUTINE findExactZeros
