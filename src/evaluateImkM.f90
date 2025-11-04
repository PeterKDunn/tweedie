SUBROUTINE evaluateImkM(i, t, f, df, m)
  ! Evaluate  Im k(t) - m*pi - pi/2  or  Im k(t) - m*pi  (CDF) for 
  ! finding the zeros of the integrand, for the PDf and CDF
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE evaluateImk(i, t_in, ImkM_out)
      ! Evaluate Im k(t)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)        :: i
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t_in
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: ImkM_out
    END SUBROUTINE evaluateImk
  END INTERFACE
  
  
  INTEGER(C_INT), INTENT(IN)        :: i
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
  INTEGER(C_INT), INTENT(IN)        :: m

  REAL(KIND=C_DOUBLE)               :: pi, Imk_val


  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  CALL evaluateImk(i, t, Imk_val)
  
  ! The expression depends on whether we are working with the PDF or the CDF.
  ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
  ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m pi/y;
  !       the CDF has integrand zeros at Im k(t) =        m pi/y.
!  WRITE(*,*) "CDF:", Cpdf
  IF (Cpdf) THEN
    f = Imk_val - REAL(m, KIND=C_DOUBLE) * pi - pi/2.0_C_DOUBLE
  ELSE
    f = Imk_val - REAL(m, KIND=C_DOUBLE) * pi
  END IF
  CALL evaluateImkd(i, t, df)
  
END SUBROUTINE evaluateImkM
