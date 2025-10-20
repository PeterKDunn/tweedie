SUBROUTINE rtsafe(i, funcd, x1, x2, xstart, xacc, root)
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  INTERFACE
      SUBROUTINE funcd_signature(i, x, f, df)
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

        IMPLICIT NONE
        REAL(KIND=C_DOUBLE), INTENT(IN)    :: x
        INTEGER(C_INT), INTENT(IN)         :: i
        REAL(KIND=C_DOUBLE), INTENT(OUT)   :: f, df

      END SUBROUTINE funcd_signature
  END INTERFACE
  
  PROCEDURE(funcd_signature) :: funcd
 

  ! Inputs
  REAL(KIND=C_DOUBLE), INTENT(IN)  :: x1, x2, xstart, xacc
  INTEGER(C_INT), INTENT(IN)       :: i

  ! Output (Function result)
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: root
  
  INTEGER, PARAMETER    :: JMAX = 50
  INTEGER               :: j
  REAL(KIND=C_DOUBLE)   :: f, df, dx, rootTMP
  REAL(KIND=C_DOUBLE)   :: rtsafeTMP, fh, fl, xh, xl
  
  ! Initialize root and boundaries
  xl = x1
  xh = x2
  rootTMP = xstart
  
  ! Get function values at boundaries
  CALL funcd(i, xl, fl, df)
  CALL funcd(i, xh, fh, df)
  
  ! --- Main Loop ---
  
  DO j = 1, JMAX    

    CALL funcd(i, rootTMP, f, df)
    
    ! Compute Newton step
    dx = f / df
    rtsafeTMP = rootTMP - dx
    
    ! Ensure step stays inside bracket
    IF (rtsafeTMP < xl .OR. rtsafeTMP > xh .OR. DABS(df) < 1.0E-12_C_DOUBLE) THEN
        rtsafeTMP = (xl + xh)/2.0_C_DOUBLE
    END IF
    
    dx = rtsafeTMP - rootTMP
    rootTMP = rtsafeTMP
    
    ! Check convergence
    IF (DABS(dx) < xacc) EXIT
    
    ! Update brackets using new f
    CALL funcd(i, rootTMP, f, df)
    IF (f < 0.0_C_DOUBLE) THEN
        xl = rootTMP
    ELSE
        xh = rootTMP
    END IF
END DO
root = rootTMP

  
END SUBROUTINE rtsafe
