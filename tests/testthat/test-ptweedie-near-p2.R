# ptweedie() for p just above 2 (2026-11). There xi = mu^(2-p)/{phi (p-2)}
# is huge whatever the dispersion (xi -> Inf as p -> 2), so the rule
# 'xi > 1e7: use the saddlepoint approximation' used it even for large phi,
# where it is poor: silent errors of 5e-4 (absolute) for phi = 0.5 and
# p = 2 + 1e-7, up to 4e-1 (relative) for phi = 5. The saddlepoint is now
# used only when phi mu^(p-2) <= 1e-5 as well; otherwise the density from
# dtweedie() is integrated.
#
# The reference: the CDF is smooth in p across p = 2, so 
#   F(y; 2 + e) = 2 F_gamma(y) - F(y; 2 - e) + O(e^2),
# and for 1 < p < 2 the series is accurate. With e = 1e-8, the O(e^2) term
# is negligible.

test_that("ptweedie() for p = 2 + 1e-8 is close to the reference", {
  e  <- 1e-8
  for (phi in c(1e-3, 0.5)) {
    y   <- qgamma(c(1e-6, 0.3, 0.7), shape = 1/phi, scale = phi)
    ref <- 2 * pgamma(y, shape = 1/phi, scale = phi) - 
           ptweedie(y, mu = 1, phi = phi, power = 2 - e)
    P   <- expect_silent( ptweedie(y, mu = 1, phi = phi, power = 2 + e) )
    # (before 2026-11: relative errors of 6e-7 for phi = 1e-3, 3e-2 for phi = 0.5)
    expect_lt( max(abs(P / ref - 1)), 
               if (phi < 0.01) 1e-12 else 2e-8 )
  }
})

test_that("ptweedie_regime() sends p near 2 with phis > 1e-5 to the integration", {
  r <- tweedie:::ptweedie_regime(c(1, 1), mu = 1, phi = c(1e-6, 0.5), power = 2 + 1e-8)
  expect_equal(r, c("saddlepoint", "integration (p near 2)"))
  # Away from p = 2, xi > 1e7 still means a small phis: the saddlepoint, as before
  r <- tweedie:::ptweedie_regime(1, mu = 1, phi = 1e-8, power = 4)
  expect_equal(r, "saddlepoint")
})

test_that("p very close to 2 with phi large: flagged, not silent", {
  # dtweedie() itself fails there, so the integration fails too, and the
  # saddlepoint approximation is used, with a warning
  expect_warning( ptweedie(1, mu = 1, phi = 5, power = 2 + 1e-9), 
                  "could not be integrated" )
})
