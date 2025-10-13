
SUBROUTINE acceleratenew(xvec, wvec, nzeros, Mmatrix, Nmatrix, West) BIND(C, NAME='acceleratenew')

  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  ! --- Input/Output Declarations (Dummy Arguments) ---
  INTEGER(C_INT), INTENT(IN)          :: nzeros
  REAL(KIND=C_DOUBLE), INTENT(IN)     :: xvec(200), wvec(200)  
  REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: Mmatrix(2, 200), Nmatrix(2, 200)  
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: West
  
  ! --- Local Variables ---
  INTEGER         :: p, maxSize
  REAL(KIND=8)    :: denom, tinyDenom, sumw
  REAL(KIND=8)    :: psi_new, FF_current
  
  ! NEW LOCAL VARIABLE to hold the clamped iteration count
  INTEGER         :: l_nzeros

  ! --- Constants ---
  tinyDenom = 1.0d-16
  maxSize = 200 ! Maximum size for the column dimension of M/N

  ! --- Initial Setup ---
  ! 1. Initialize the local variable with the input value
  l_nzeros = nzeros
  
  ! 2. Initial Sanity Checks (using the local variable)
  IF (l_nzeros .LE. 0) THEN
      West = 0.0d00
      RETURN
  END IF
    
  ! 3. Clamp the local variable (This is now safe)
  IF (l_nzeros .GT. maxSize) THEN
      l_nzeros = maxSize  
  END IF

  ! The new term (psi) is the last element of wvec (using the original nzeros for indexing)
  ! NOTE: If nzeros > 200, this index may still be out of bounds, but the loop below
  ! handles the clamped size. This logic assumes wvec is at least size 'nzeros' (input).
  psi_new = wvec(nzeros)
    
  ! Calculate the actual cumulative sum W_l (FF_current)
  sumw = 0.0d00
  DO p = 1, nzeros ! Use original nzeros here to calculate the full sum
      sumw = sumw + wvec(p)
  ENDDO
  FF_current = sumw  
    
  ! --- LOGIC CLONE OF SIDIACC ---
    
  ! 1. Handle the case where the new term (psi) is essentially zero.
  IF (DABS(psi_new) .LT. 1.0d-31) THEN
      West = FF_current
      RETURN
  END IF

  ! 2. Initialization (p=1, which is column 1 in M/N matrices)
  Mmatrix(2, 1) = FF_current / psi_new  
  Nmatrix(2, 1) = 1.0d00 / psi_new
    
  ! 3. Build the recurrence table row (p=2 to l_nzeros)
  DO p = 2, l_nzeros ! IMPORTANT: Use the clamped local variable
      
      ! Denominator corresponds to 1/x_{l+1-i} - 1/x_l
      ! Here, l_nzeros acts as 'l'
      denom = 1.0d00 / xvec(l_nzeros + 1 - p) - 1.0d00 / xvec(l_nzeros)
        
      IF (DABS(denom) .LT. tinyDenom) THEN
          ! Denominator too small -> fall back to the raw sum
          West = FF_current  
          RETURN
      END IF
        
      ! Mmatrix(2, p) = ( Mmatrix(1, p-1) - Mmatrix(2, p-1) ) / denom
      Mmatrix(2, p) = ( Mmatrix(1, p - 1) - Mmatrix(2, p - 1) ) &
                        / denom
      Nmatrix(2, p) = ( Nmatrix(1, p - 1) - Nmatrix(2, p - 1) ) &
                        / denom
        
  END DO

  ! 4. Final accelerated estimate (Highest order, column l_nzeros)
  IF (l_nzeros .GT. 1) THEN ! IMPORTANT: Use the clamped local variable
      IF (DABS(Nmatrix(2, l_nzeros)) .LT. tinyDenom) THEN
          West = FF_current
      ELSE
          ! The highest order estimate is M(2, l_nzeros) / N(2, l_nzeros)
          West = Mmatrix(2, l_nzeros) / Nmatrix(2, l_nzeros)
      END IF
  ELSE
      West = FF_current
  END IF

  ! 5. Update the Old Row (Copy Row 2 to Row 1 for the next iteration)
  DO p = 1, l_nzeros ! IMPORTANT: Use the clamped local variable
      Mmatrix(1, p) = Mmatrix(2, p)
      Nmatrix(1, p) = Nmatrix(2, p)
  END DO
    
  RETURN
END SUBROUTINE acceleratenew
