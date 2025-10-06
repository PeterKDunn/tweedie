      SUBROUTINE accelerate(xvec, wvec, nzeros, pmax, West)
*
*     SIDI_ACCEL (simple, robust implementation)
*     xvec(1..nzeros+1)  zeros: x0,x1,...,x_n
*     wvec(1..nzeros)    areas: integral_{x_{s}}^{x_{s+1}} f(t) dt
*     nzeros >= 1
*     pmax   requested acceleration order (non-negative)
*
      IMPLICIT NONE

      INTEGER nzeros, pmax
      DOUBLE PRECISION xvec(400), wvec(400), West
      DOUBLE PRECISION Mmatrix(400, 50), Nmatrix(400, 50)
      DOUBLE PRECISION denom, sumw, eps, tinyDenom
      DOUBLE PRECISION sumw_temp
      INTEGER s, p, pmax_use, maxSize

*     M and N are indexed as M(s, p)

*     Tolerances
      eps = 1.0d-100
      tinyDenom = 1.0d-16
      maxSize = 400

*     Default
      West = 0.0d00

*     Sanity checks
      IF (nzeros .LE. 0) THEN
         West = 0.0d00
         RETURN
      ENDIF
      
*     Limit nzeros to array size (maxSize = 200)
      IF (nzeros .GT. maxSize - 1) THEN
          nzeros = maxSize - 1
      ENDIF

      IF (nzeros .EQ. 1) THEN
         West = wvec(1)
         RETURN
      ENDIF

*     Clamp pmax to available data and matrix width (second dim 50)
      IF (pmax .GT. nzeros - 1) THEN
         pmax_use = nzeros - 1
      ELSE
         pmax_use = pmax
      ENDIF
      IF (pmax_use .GT. 48) pmax_use = 48   ! keep p+1 <= 49 inside 50 cols

*     Unaccelerated fallback sum
      sumw = 0.0d00
      DO s = 1, nzeros
         sumw = sumw + wvec(s)
      ENDDO

*     Initialize M(:,1) and N(:,1) (Partial Sums - W(s, 0) = M(s,1)/N(s,1))
      sumw_temp = 0.0d00
      DO s = 1, nzeros
          sumw_temp = sumw_temp + wvec(s)
          Mmatrix(s,1) = sumw_temp
          Nmatrix(s,1) = 1.0d00
      ENDDO

*     If any w is essentially zero, Sidi acceleration is unreliable:
      DO s = 1, nzeros
         IF (DABS(wvec(s)) .LT. eps) THEN
            West = sumw
            RETURN
         ENDIF
      ENDDO

*     Build recurrence:
*       for p = 1..pmax_use
*         for s = 1..(nzeros - p)
*           denom = 1/x_s - 1/x_{s+p}   (xvec must be zeros in increasing order)
*           M(s,p+1) = ( M(s,p) - M(s+1,p) ) / denom
*           N(s,p+1) = ( N(s,p) - N(s+1,p) ) / denom
      DO p = 1, pmax_use
         DO s = 1, nzeros - p
            denom = 1.0d00 / xvec(s) - 1.0d00 / xvec(s + p)
            IF (DABS(denom) .LT. tinyDenom) THEN
*              denominator too small -> fallback
               West = sumw
               RETURN
            ENDIF
            Mmatrix(s, p+1) = ( Mmatrix(s, p) - Mmatrix(s+1, p) ) 
     &                         / denom
            Nmatrix(s, p+1) = ( Nmatrix(s, p) - Nmatrix(s+1, p) ) 
     &                         / denom
         ENDDO
      ENDDO

*     Final check and compute accelerated estimate (col index = pmax_use+1)
      IF (DABS(Nmatrix(1, pmax_use + 1)) .LT. tinyDenom) THEN
         West = sumw
      ELSE
         West = Mmatrix(1, pmax_use + 1) / Nmatrix(1, pmax_use + 1)
      ENDIF

      RETURN
      END
