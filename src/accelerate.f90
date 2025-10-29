SUBROUTINE accelerate(xvec, wvec, nzeros, Mmatrix, Nmatrix, West)
  ! Accelerate the convergence of the infinite integral
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  IMPLICIT NONE

  INTEGER, INTENT(IN)                 :: nzeros
  REAL(KIND=C_DOUBLE), INTENT(IN)     :: xvec(200), wvec(200)
  REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: Mmatrix(2,200), Nmatrix(2,200)
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: West

  ! --- Local variables ---
  INTEGER                         :: p, l_nzeros, maxSize
  REAL(KIND=C_DOUBLE)             :: denom, psi_new, FF_current, sumw
  REAL(KIND=C_DOUBLE)             :: tinyDenom, scale_denom
  REAL(KIND=C_DOUBLE)             :: inv_x_i, inv_x_l
  REAL(KIND=C_DOUBLE)             :: s, maxinv, invx
  REAL(KIND=C_DOUBLE), PARAMETER  :: SAFETY_SCALE = 1.0D-12
  REAL(KIND=C_DOUBLE), PARAMETER  :: HUGE_LIMIT   = 1.0D300
  REAL(KIND=C_DOUBLE)             :: xscaled(200)

  ! --- Constants; initialization ---
  maxSize  = 200
  l_nzeros = MIN(nzeros, maxSize)

  IF (l_nzeros .LE. 0) THEN
     West = 0.0D0
     RETURN
  END IF

  ! --- Rescaling, to improve numerical conditioning ---
  maxinv = 0.0D0
  DO p = 1, l_nzeros
     IF (xvec(p) .EQ. 0.0E0_C_DOUBLE) THEN
        West = 0.0E0_C_DOUBLE
        RETURN
     END IF
     invx = DABS(1.0E0_C_DOUBLE / xvec(p))
     IF (invx .GT. maxinv) maxinv = invx
  END DO
  s = MAX(1.0E0_C_DOUBLE, 1.0E0_C_DOUBLE / maxinv)

  DO p = 1, l_nzeros
     xscaled(p) = xvec(p) / s
  END DO


  ! --- BEGIN: Main algorithm ---
  psi_new = wvec(nzeros)

  sumw = 0.0E0_C_DOUBLE
  DO p = 1, nzeros
     sumw = sumw + wvec(p)
  END DO
  FF_current = sumw

  IF (DABS(psi_new) .LT. 1.0E-31_C_DOUBLE) THEN
     West = FF_current
     RETURN
  END IF

  Mmatrix(2, 1) = FF_current / psi_new
  Nmatrix(2, 1) = 1.0E0_C_DOUBLE / psi_new

  DO p = 2, l_nzeros
     inv_x_l = 1.0E0_C_DOUBLE / xscaled(l_nzeros)
     inv_x_i = 1.0E0_C_DOUBLE / xscaled(l_nzeros + 1 - p)
     denom   = inv_x_i - inv_x_l
     scale_denom = MAX( DABS(inv_x_i),  & 
                        DABS(inv_x_l),  &
                        1.0E0_C_DOUBLE )
     tinyDenom   = SAFETY_SCALE * scale_denom

     IF ( (.NOT.(denom .EQ. denom)) .OR.  & 
          (DABS(denom) .LT. tinyDenom) ) THEN
        West = FF_current
        RETURN
     END IF

     Mmatrix(2, p) = (Mmatrix(1, p - 1) - Mmatrix(2, p - 1)) / denom
     Nmatrix(2, p) = (Nmatrix(1, p - 1) - Nmatrix(2, p - 1)) / denom

     IF (.NOT.(Mmatrix(2, p) .EQ. Mmatrix(2, p)) .OR.   &
         .NOT.(Nmatrix(2, p) .EQ. Nmatrix(2, p))) THEN
        West = FF_current
        RETURN
     END IF

     IF ( (DABS(Mmatrix(2, p)) .GT. HUGE_LIMIT) .OR.    &
          (DABS(Nmatrix(2, p)) .GT. HUGE_LIMIT) ) THEN
        West = FF_current
        RETURN
     END IF
  END DO
  ! --- END: Main algorithm ---


  ! Final accelerated estimate
  IF (l_nzeros > 1) THEN
    IF (.NOT.(Nmatrix(2, l_nzeros) .EQ. Nmatrix(2, l_nzeros))) THEN
       ! NaN detected, use previous column
       West = Mmatrix(2, l_nzeros - 1) / Nmatrix(2, l_nzeros - 1)
       RETURN
    END IF
    IF (DABS(Nmatrix(2, l_nzeros)) .LT. tinyDenom .OR. &
        DABS(Mmatrix(2, l_nzeros)) .GT. 1.0E300_C_DOUBLE) THEN
       ! Unstable last column: use lower order
       West = Mmatrix(2, l_nzeros-1) / Nmatrix(2, l_nzeros - 1)
       RETURN
    ELSE
       West = Mmatrix(2, l_nzeros) / Nmatrix(2, l_nzeros)
    END IF
  ELSE
    West = FF_current
  END IF


  ! --- Copy row 2 to row 1 ---
  DO p = 1, l_nzeros
     Mmatrix(1, p) = Mmatrix(2, p)
     Nmatrix(1, p) = Nmatrix(2, p)
  END DO

  RETURN
END SUBROUTINE accelerate
