# For p > 2 and very small dispersion (xi = mu^(2-p) / {phi (p-2)} > 1e5),
# ptweedie() uses the Lugannani-Rice saddlepoint approximation
# (ptweedie_saddle), since 3.2.1.

test_that("ptweedie_saddle() matches the exact inverse Gaussian (p = 3) at small dispersion", {
  for (xi in c(1e6, 1e8, 1e10)) {
    phi <- 1 / xi                            # p = 3, mu = 1
    sd  <- sqrt(phi)
    y   <- 1 + c(-8, -4, -1, 0, 1, 4, 8, 20) * sd
    for (lt in c(TRUE, FALSE)) {
      ref <- statmod::pinvgauss(y, mean = 1, dispersion = phi, lower.tail = lt, log.p = TRUE)
      lr  <- tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = 3,
                                       lower.tail = lt, log.p = TRUE)
      # compare probabilities on the relative scale
      expect_lt(max(abs(expm1(lr - ref))), 1e-8)
    }
  }
})

test_that("ptweedie_saddle() is finite, increasing and smooth near y = mu", {
  # Near y = mu a Taylor series is used for 1/w - 1/u (|log(y/mu)| < 1e-5).
  # Check there is no jump where it switches: across a step h the change in F
  # must be close to (density) x h.
  for (p in c(2.05, 3, 8)) {
    phi <- 1 / (1e7 * (p - 2))
    y   <- exp(c(-1e-5 - 1e-9, -1e-5 + 1e-9, -1e-12, 0, 1e-12, 1e-5 - 1e-9, 1e-5 + 1e-9))
    F   <- tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p)
    expect_true(all(is.finite(F)))
    expect_true(all(diff(F) > 0))
    expect_equal(F[4], 0.5, tolerance = 1e-3)
    for (i in c(1, 6)) {
      slope <- (F[i + 1] - F[i]) / (y[i + 1] - y[i])
      dens  <- dtweedie_saddle(y = (y[i] + y[i + 1]) / 2, mu = 1, phi = phi, power = p)
      expect_equal(slope, as.numeric(dens), tolerance = 1e-3)
    }
  }
})

test_that("ptweedie() uses the saddlepoint approximation for xi > 1e5 (p > 2)", {
  p   <- 2.5
  xi  <- 1e7
  phi <- 1 / (xi * (p - 2))
  y   <- 1 + c(-6, -2, 0, 2, 6) * sqrt(phi)
  expect_silent(F <- ptweedie(y, mu = 1, phi = phi, power = p))
  expect_equal(F, tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p))
  expect_silent(S <- ptweedie(y, mu = 1, phi = phi, power = p, lower.tail = FALSE, log.p = TRUE))
  expect_equal(S, tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p,
                                            lower.tail = FALSE, log.p = TRUE))
})

test_that("ptweedie() for p > 2 is close to continuous at the switch, xi = 1e5", {
  # At xi = 1e5 the inversion (with tail integration) is still used; just
  # above it, the saddlepoint approximation, whose error there is about 1e-7.
  for (p in c(2.05, 4)) {
    phi <- 1 / (1e5 * (p - 2))
    y   <- 1 + c(-4, -1, 0, 1, 4) * sqrt(phi)
    inv <- ptweedie(y, mu = 1, phi = phi, power = p)
    sp  <- tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p)
    expect_lt(max(abs(sp / inv - 1)), 1e-6)
  }
})
