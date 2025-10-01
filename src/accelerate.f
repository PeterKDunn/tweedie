      SUBROUTINE accelerate(xvec, wvec, nzeros, pmax, West)
*
*     SIDI_ACCEL: Sidi extrapolation/acceleration for oscillatory integrals
*     using exact zeros.
*
*     IN:
*        x(nzeros+1) : DOUBLE PRECISION array of zeros of f(t)
*        w(nzeros)   : DOUBLE PRECISION array of integrals over [x_s, x_{s+1}]
*        nzeros      : INTEGER, number of intervals
*        pmax        : INTEGER, maximum acceleration order
*
*     OUT:
*        W           : DOUBLE PRECISION, accelerated integral estimate
*
      IMPLICIT NONE
      DOUBLE PRECISION xvec(200), wvec(200), West
      DOUBLE PRECISION Mmatrix(200, 50)
      DOUBLE PRECISION Nmatrix(200, 50)
* the second dimension should be accMax + 2
      DOUBLE PRECISION denom
      INTEGER nzeros, pmax
      INTEGER s, p


*     ERROR checks: If the latest region is zero, just return; we
*                   won't be able to improve on the estimate, West
      IF (wvec(1) .EQ. 0.0d00) THEN
        West = 0.0d0
        write(*,*) "FIRST interval has zero area!"
        RETURN
      ENDIF
      IF (wvec(2) .EQ. 0.0d00) THEN
        write(*,*) "SECOND interval has zero area!"
*       So just return West
        RETURN
      ENDIF
      
      IF (wvec(nzeros) .EQ. 0.0d00) THEN
        West = Mmatrix(1, 1)  
        write(*,*) "LATEST interval has zero area!"
*       Return the first interval
         RETURN
      ENDIF

*     Trivial case: only one interval
      IF (nzeros .LE. 1) THEN
         West = wvec(1)
         write(*,*) "Only one interval:"
         write(*,*) "returning unaccelerated estimate."
         RETURN
      ENDIF

*     Initialize M and N arrays
*     Note: Index 1 corresponds to paper's "order -1"
      DO s = 1, nzeros
         Mmatrix(s, 1) = wvec(s)
*        Corresponds to M_{-1}^{(s-1)}
         Nmatrix(s, 1) = 1.0d00 / wvec(s)
      END DO

*     Recursive acceleration
      IF (pmax > nzeros - 1) pmax = nzeros - 1
      DO p = 1, (pmax + 1)
         DO s = 1, nzeros - p
            denom = 1.0d00 / xvec(s) - 1.0d00 / xvec(s + p)
            Mmatrix(s, p + 1) = (Mmatrix(s, p) - 
     &                           Mmatrix(s + 1, p)) / denom
            Nmatrix(s, p + 1) = (Nmatrix(s, p) - 
     &                           Nmatrix(s + 1, p)) / denom
         ENDDO
      ENDDO

*     Compute accelerated integral estimate
      West = Mmatrix(1, pmax + 1) / Nmatrix(1, pmax + 1)

      RETURN
      END
