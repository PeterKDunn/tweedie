# Cases that previously never returned (infinite loop in improveKZeroBounds,
# caused by mmax being computed too large for the PDF). If these regress, the
# test will hang rather than fail.

test_that("PDF with kmax > pi (p = 2.05) matches a high-precision reference", {
  # (This case hung when called on the unscaled scale; via the public API it
  # is an accuracy check on the mmax change.)
  # Reference: 40-digit mpmath integration of the inversion integral
  expect_equal(dtweedie_inversion(0.005, mu = 0.01, phi = 0.01, power = 2.05),
               1.85832121185e-08, tolerance = 1e-9)
})

test_that("PDF with kmax > pi returns (p = 1.3), matching the series method", {
  expect_equal(dtweedie_inversion(1e-5, mu = 0.01, phi = 0.01, power = 1.3),
               as.numeric(dtweedie_series(1e-5, power = 1.3, mu = 0.01, phi = 0.01)),
               tolerance = 1e-9)
})

test_that("PDF far in the upper tail (p > 2) returns ~0 instead of hanging", {
  # Saddlepoint log-densities here are below -3000, so the density underflows
  expect_equal(dtweedie_inversion(0.1, mu = 0.01, phi = 1,  power = 3.5), 0)
  expect_equal(dtweedie_inversion(0.1, mu = 0.01, phi = 1,  power = 5),   0)
  expect_equal(dtweedie_inversion(0.1, mu = 0.01, phi = 10, power = 5),   0)
})
