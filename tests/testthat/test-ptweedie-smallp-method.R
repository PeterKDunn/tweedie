# For 1 < p < 2, ptweedie() uses the series when lambda <= 1e5 (exact: it sums
# positive terms), and the inversion otherwise.

test_that("ptweedie: large point mass, p near 1 (inversion was silently wrong)", {
  # lambda ~ 0.094; the inversion gave 0.91893 here, with no warning
  y <- c(0.001, 0.01, 0.5)
  s <- ptweedie_series(y, mu = 10, phi = 100, power = 1.05)
  expect_silent(v <- ptweedie(y, mu = 10, phi = 100, power = 1.05))
  expect_equal(v, s, tolerance = 1e-14)
  expect_lt(abs(v[1] - 0.9104504), 1e-6)
})

test_that("ptweedie: tiny q for p near 2 (Lacerda, 2017)", {
  # Reference: 50-digit sum of the series (mpmath). Not just Pr(Y = 0): the
  # gamma shapes are tiny (about 0.02 j), so the continuous part contributes.
  v <- ptweedie(7.709933e-308, mu = 1.017691e+01, phi = 4.55, power = 1.98)
  expect_equal(v, 1.0019896546587636e-5, tolerance = 1e-12)
})

test_that("ptweedie: very large lambda still uses the inversion, without warning", {
  # lambda = 1e6: beyond the series' range; the inversion is reliable here
  expect_silent(v <- ptweedie(1, mu = 1, phi = 1e-4, power = 1.99))
  expect_equal(v, ptweedie_inversion(1, mu = 1, phi = 1e-4, power = 1.99))
})

test_that("ptweedie: agrees with pchisq for p = 1.5 in both tails and on the log scale", {
  y <- seq(0.05, 15, length = 40)
  expect_equal(ptweedie(y, mu = 1, phi = 4, power = 1.5), pchisq(y, df = 0, ncp = 1), tolerance = 1e-14)
  expect_equal(ptweedie(y, mu = 1, phi = 4, power = 1.5, log.p = TRUE),
               pchisq(y, df = 0, ncp = 1, log.p = TRUE), tolerance = 1e-13)
  expect_equal(ptweedie(y[1:20], mu = 1, phi = 4, power = 1.5, lower.tail = FALSE),
               pchisq(y[1:20], df = 0, ncp = 1, lower.tail = FALSE), tolerance = 1e-12)
})
