# For 1 < p < 2, ptweedie() uses the series at every lambda (since 3.2.1).
# Before, it used the inversion for lambda > 1e9, which for lambda >= 1e10
# returned about 0.5 (with a warning) for y below the mean.

test_that("ptweedie() is accurate for very large lambda (1 < p < 2)", {
  p   <- 1.5
  lam <- 1e10
  phi <- 1 / (lam * (2 - p))
  sd  <- sqrt(phi)                       # mu = 1
  y   <- 1 + c(-2, -1, 1, 2) * sd
  ## Reference: the Lugannani-Rice approximation, whose error is far below
  ## the tolerance here (lambda = 1e10), computed without cancellation
  e  <- log(y)
  th <- -expm1((1 - p) * e) / (phi * (p - 1))
  d  <- 2 * (expm1((2 - p) * e) - (2 - p) * expm1(e)) / ((1 - p) * (2 - p))
  w  <- sign(y - 1) * sqrt(d / phi)
  u  <- th * sqrt(phi * y^p)
  ref <- pnorm(w) + dnorm(w) * (1 / w - 1 / u)
  expect_silent(F <- ptweedie(y, mu = 1, phi = phi, power = p))
  expect_equal(F, ref, tolerance = 1e-8)
})
