# The inversion computes F = 1/2 -/+ I/pi, so its error is ABSOLUTE: small
# tail probabilities have poor relative accuracy. ptweedie_inversion() flags
# them (small_tail, with details = TRUE); ptweedie() recomputes them by
# integrating the density. Values are checked ONE BY ONE (not with
# expect_equal() on a whole vector, whose mean relative difference hides a
# poor tail value among accurate ones).
#
# Reference: the inverse Gaussian (p = 3), with IGexact = FALSE so that the
# inversion is actually used.

mu  <- 1.4
phi <- 0.74

test_that("inversion: lower tail accurate in the body; small tails flagged", {
  y <- c(0.03, 0.05, 0.1, 0.5, 1, 2, 5)
  F <- statmod::pinvgauss(y, mean = mu, dispersion = phi)
  out <- ptweedie_inversion(y, mu = mu, phi = phi, power = 3,
                            IGexact = FALSE, details = TRUE)

  # Flagged exactly where the returned value is below 1e-8 (here y = 0.03)
  expect_identical(out$small_tail, out$cdf < 1e-8)
  expect_identical(out$small_tail, y == 0.03)

  # Not flagged: good RELATIVE accuracy, value by value
  ok <- !out$small_tail
  expect_true(all(abs(out$cdf[ok] / F[ok] - 1) < 1e-8))
  body <- F > 1e-5
  expect_true(all(abs(out$cdf[body] / F[body] - 1) < 1e-12))

  # Flagged: still good ABSOLUTE accuracy
  expect_true(all(abs(out$cdf[!ok] - F[!ok]) < 1e-14))
})

test_that("inversion: upper tail flags small tail probabilities", {
  y <- c(5, 20, 40, 60)
  S <- statmod::pinvgauss(y, mean = mu, dispersion = phi, lower.tail = FALSE)
  out <- ptweedie_inversion(y, mu = mu, phi = phi, power = 3, IGexact = FALSE,
                            lower.tail = FALSE, details = TRUE)
  expect_identical(out$small_tail, out$cdf < 1e-8)
  expect_identical(out$small_tail, y == 60)
  expect_true(all(abs(out$cdf - S) < 1e-14))
  body <- S > 1e-5
  expect_true(all(abs(out$cdf[body] / S[body] - 1) < 1e-12))
})

test_that("inversion: the flag is set on the probability scale when log.p = TRUE", {
  y   <- c(0.03, 0.05, 1)
  out <- ptweedie_inversion(y, mu = mu, phi = phi, power = 3, IGexact = FALSE,
                            log.p = TRUE, details = TRUE)
  expect_identical(out$small_tail, c(TRUE, FALSE, FALSE))
})

test_that("inversion: exact special cases are never flagged; no warning is given", {
  out <- ptweedie_inversion(c(0.03, 60), mu = mu, phi = phi, power = 3,
                            details = TRUE)                  # IGexact = TRUE
  expect_identical(out$small_tail, c(FALSE, FALSE))
  expect_no_warning(
    ptweedie_inversion(c(0.03, 0.05), mu = mu, phi = phi, power = 3, IGexact = FALSE)
  )
})

test_that("the tail repair used by ptweedie() has relative accuracy, value by value", {
  for (q in c(0.03, 0.05)) {
    ti <- tweedie:::ptweedie_tail_integrate(q, mu, phi, power = 3, IGexact = FALSE)
    lF <- statmod::pinvgauss(q, mean = mu, dispersion = phi, log.p = TRUE)
    expect_true(ti$ok)
    expect_lt(abs(exp(ti$logp - lF) - 1), 1e-12)
  }
  for (q in c(40, 60)) {
    ti <- tweedie:::ptweedie_tail_integrate(q, mu, phi, power = 3,
                                            lower.tail = FALSE, IGexact = FALSE)
    lS <- statmod::pinvgauss(q, mean = mu, dispersion = phi,
                             lower.tail = FALSE, log.p = TRUE)
    expect_true(ti$ok)
    expect_lt(abs(exp(ti$logp - lS) - 1), 1e-12)
  }
})
