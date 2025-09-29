      SUBROUTINE findRek(p, mu, phi, y)

      DOUBLE PRFECISION lambda, p, mu, phi, y
      DOUBLE PRFECISION kmax, tmax, startPoint
      DOUBLE PRFECISION front, alpha, omega, t
      INTEGER mmax

      pi = dacos(-1.0d00)

      write(*,*) "- Finding Re(k(t))", y

      front = mu ** (2.0d00 - p) / ( phi*(2.0d00 - p))
      omega = datan( (1.0d00 - p) * t * phi)/
     &               (mu ** (1.0d00 - p) )
      alpha = (2.0d00 - p)/(1.0d00 - p)

      Imk = front *
     &      ( cos(omega * alpha)/(cos(omega)**alpha) - 1.0d00)

      RETURN
      END