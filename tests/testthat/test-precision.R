# Cases where Re k(t), Im k(t) were evaluated inaccurately for large
# |tanArg| = c*t (DATAN rounding to exactly -pi/2, so DCOS(omega) was ~6e-17
# instead of ~1/|tanArg|), and/or the first integration region spanned many
# decades beyond t ~ 1/c and was under-resolved by a single Gauss rule.

test_that("p > 2 densities that were negative now match the series method", {
  cases <- list(c(500, 50, 10, 5), c(500, 50, 1, 5), c(75, 50, 10, 8),
                c(25, 50, 10, 8),  c(45, 50, 10, 8))
  for (z in cases) {
    inv <- dtweedie_inversion(z[1], mu = z[2], phi = z[3], power = z[4])
    ser <- as.numeric(dtweedie_series(z[1], power = z[4], mu = z[2], phi = z[3]))
    expect_gt(inv, 0)
    expect_equal(inv, ser, tolerance = 1e-9)
  }
})

test_that("extreme lower tail (p = 8) is ~0 and reported as converged", {
  # Chernoff bound: P(Y <= 0.05) < exp(-152000), so the CDF is 0 in double precision
  out <- ptweedie_inversion(0.05, mu = 50, phi = 10, power = 8, details = TRUE)
  expect_equal(out$cdf, 0, tolerance = 1e-12)
  expect_equal(out$exitstatus, 0L)
})

test_that("upper-tail CDFs for large p match an independent reference", {
  # Reference: independent double-precision Gauss-Legendre integration of the
  # inversion integral (graded panels near 0, panels of 1/8 period beyond)
  expect_equal(ptweedie_inversion(500, mu = 50, phi = 10, power = 8), 0.99952702, tolerance = 1e-7)
  expect_equal(ptweedie_inversion(500, mu = 50, phi = 1,  power = 8), 0.99934221, tolerance = 1e-7)
  expect_equal(ptweedie_inversion(75,  mu = 50, phi = 1,  power = 8), 0.99661162, tolerance = 1e-7)
})
