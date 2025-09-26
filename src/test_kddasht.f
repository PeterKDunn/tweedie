      program test_kddasht
      double precision p, phi, mu, x, ddrl, ddim

      ! Example values from your log that caused segfault
      p = 3.0d0
      phi = 1.0d0
      mu = 0.1d0
      x = 1.0d0

      call kddasht(p, phi, mu, x, ddrl, ddim)

      print *, "ddrl =", ddrl
      print *, "ddim =", ddim

      end
