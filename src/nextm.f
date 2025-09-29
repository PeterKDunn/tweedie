
      SUBROUTINE  nextm( tmax, mmax, zero2, m, firstm,
     &                   zlo, zhi, zero )
*     This subroutine advances to the next value of  m  in the
*     case  p > 2, and  y < mu  and  kmax > pi.  It also returns,
*     when appropriate, adjusted bounds on the next zero.

*     IN:    tmax, mmax, zero2, m, firstm, zlo
*     OUT:   m, firstm, zlo, zhi, zero

      double precision  tmax, zero2, zlo, zhi, zero
      integer  mmax, m, firstm

*     Can get tricky, as need to be careful when  Im(k)  turns back down.
      if ( m .LT. mmax ) then
         if ( firstm .EQ. 1 ) then
            m = m + 1
            zhi = tmax
         else
            m = m - 1
            zlo = max( zlo, tmax )
         endif
      elseif ( m .EQ. mmax ) then
         if ( firstm .EQ. 1 ) then
            firstm = firstm + 1
            zero = tmax + (tmax - zero2)
            zlo = tmax
         else
            m = m - 1
*            zlo = tmax
         endif
      endif

      return
      end
