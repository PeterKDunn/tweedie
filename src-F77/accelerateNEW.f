      SUBROUTINE accelerateNEW(xvec, wvec, nzeros, 
     &                         Mmatrix, Nmatrix, West)
*
* ACCELERATE2 - Incremental Sidi D-Transformation (SIDIACC Clone)
*
* This subroutine mimics the behavior of SIDIACC: it performs
* one step of acceleration, using the previous highest-order
* estimates (stored in Mmatrix(1,:) and Nmatrix(1,:)) and the
* newest term (wvec(nzeros)) to calculate the new set of estimates
* (stored in Mmatrix(2,:) and Nmatrix(2,:)).
*
* Inputs:
* xvec(1..nzeros) : The x-values used up to the current point.
* wvec(1..nzeros) : The individual integral areas (psi_i).
* nzeros          : The current iteration count / number of terms (l).
* * Input/Output (Persistent Table):
* Mmatrix(2, 200) : The numerator matrix (row 1 is OLD, row 2 is NEW).
* Nmatrix(2, 200) : The denominator matrix (row 1 is OLD, row 2 is NEW).
*
* Output:
* West            : The newest accelerated estimate of the infinite sum (W).
*
      IMPLICIT NONE

      ! --- Input/Output Declarations (Matching SIDIACC usage) ---
      INTEGER nzeros
      DOUBLE PRECISION xvec(*), wvec(*) 
      ! M and N arrays are assumed to be sized M(2, N_MAX), N(2, N_MAX)
      DOUBLE PRECISION Mmatrix(2, 200), Nmatrix(2, 200) 
      DOUBLE PRECISION West

      ! --- Local Variables (Matching ACCELERATE style) ---
      INTEGER p, maxSize
      DOUBLE PRECISION denom, tinyDenom, sumw
      DOUBLE PRECISION psi_new, FF_current

      ! --- Constants ---
      tinyDenom = 1.0d-16
      maxSize = 200 ! Maximum size for the column dimension of M/N

      ! --- Initial Sanity Checks ---
      IF (nzeros .LE. 0) THEN
          West = 0.0d00
          RETURN
      ENDIF
      
      ! Clamp nzeros to matrix column size limit
      IF (nzeros .GT. maxSize) THEN
          ! Since nzeros must be used to index M/N, this protects the array bounds.
          nzeros = maxSize 
      ENDIF

      ! The new term (psi) is the last element of wvec
      psi_new = wvec(nzeros)
      
      ! FF (The cumulative sum up to x_{l-1}) needs to be calculated
      sumw = 0.0d00
      DO p = 1, nzeros
          sumw = sumw + wvec(p)
      ENDDO
      FF_current = sumw - psi_new ! FF in sidiacc is the sum *up to* the previous step. 
                                  ! But sidiacc's setup is FF (sum up to x_{l})

      ! Recalculate FF_current as the actual cumulative sum W_l
      FF_current = sumw 
      
      ! --- LOGIC CLONE OF SIDIACC ---
      
      ! 1. Handle the case where the new term (psi) is essentially zero.
      IF (DABS(psi_new) .LT. 1.0d-31) THEN
          West = FF_current
          RETURN
      ENDIF

      ! 2. Initialization (p=1, which is column 1 in M/N matrices)
      ! Note: sidiacc uses M(2,1) and N(2,1) for the *newest* information.
      ! This is the non-standard initialization that makes sidiacc unique.
      Mmatrix(2, 1) = FF_current / psi_new 
      Nmatrix(2, 1) = 1.0d00 / psi_new
      
      ! 3. Build the recurrence table row (p=2 to nzeros)
      ! The index p corresponds to 'i' in the sidiacc loop.
      ! The loop calculates estimates of order p-1 up to nzeros-1.
      DO p = 2, nzeros
          
          ! Denominator corresponds to 1/x_{l+1-i} - 1/x_l
          ! Where l = nzeros and i = p. The index in xvec starts at 1.
          ! xvec(nzeros + 1 - p) is x_{l-p+1}
          ! xvec(nzeros) is x_l
          denom = 1.0d00 / xvec(nzeros + 1 - p) - 1.0d00 / xvec(nzeros)
          
          IF (DABS(denom) .LT. tinyDenom) THEN
              ! Denominator too small -> fall back to the raw sum
              West = FF_current 
              RETURN
          ENDIF
          
          ! Mmatrix(2, p) = ( Mmatrix(1, p-1) - Mmatrix(2, p-1) ) / denom
          ! (New estimate of order p-1) = (Old order p-2 - New order p-2) / denom
          Mmatrix(2, p) = ( Mmatrix(1, p - 1) - Mmatrix(2, p - 1) )
     &                    / denom
          Nmatrix(2, p) = ( Nmatrix(1, p - 1) - Nmatrix(2, p - 1) ) 
     &                    / denom
          
      ENDDO

      ! 4. Final accelerated estimate (Highest order, column nzeros)
      ! Only if nzeros > 1, otherwise it falls back to the FF_current/sumw
      IF (nzeros .GT. 1) THEN
          IF (DABS(Nmatrix(2, nzeros)) .LT. tinyDenom) THEN
              West = FF_current
          ELSE
              ! The highest order estimate is M(2, nzeros) / N(2, nzeros)
              West = Mmatrix(2, nzeros) / Nmatrix(2, nzeros)
          ENDIF
      ELSE
          West = FF_current
      ENDIF

      ! 5. Update the Old Row (Copy Row 2 to Row 1 for the next iteration)
      ! This 'drops' the first row by promoting the new row to the old row.
      DO p = 1, nzeros
          Mmatrix(1, p) = Mmatrix(2, p)
          Nmatrix(1, p) = Nmatrix(2, p)
      ENDDO
      
      RETURN
      END 
