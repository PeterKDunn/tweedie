# For p > 2 and very small dispersion, ptweedie() uses the Lugannani-Rice
# saddlepoint approximation (ptweedie_saddle), since 3.2.0: for every y when
# xi = mu^(2-p) / {phi (p-2)} > 1e7, and in the body when 1e6 < xi <= 1e7 and
# phi mu^(p-2) <= 1e-5.

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

test_that("ptweedie() uses the saddlepoint approximation for xi > 1e7 (p > 2)", {
  p   <- 2.5
  xi  <- 1e8
  phi <- 1 / (xi * (p - 2))
  y   <- 1 + c(-6, -2, 0, 2, 6) * sqrt(phi)
  expect_silent(F <- ptweedie(y, mu = 1, phi = phi, power = p))
  expect_equal(F, tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p))
  expect_silent(S <- ptweedie(y, mu = 1, phi = phi, power = p, lower.tail = FALSE, log.p = TRUE))
  expect_equal(S, tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p,
                                            lower.tail = FALSE, log.p = TRUE))
})

test_that("ptweedie() for p > 2 is close to continuous at the switch, xi = 1e7", {
  # At xi = 1e7 the inversion and integration are still used; just above it,
  # the saddlepoint approximation. Both are accurate to about 1e-9 here.
  for (p in c(2.05, 4)) {
    phi <- 1 / (1e7 * (p - 2))
    y   <- 1 + c(-8, -4, -1, 0, 1, 4, 8) * sqrt(phi)
    for (lt in c(TRUE, FALSE)) {
      a  <- ptweedie(y, mu = 1, phi = phi, power = p, lower.tail = lt)
      sp <- tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p, lower.tail = lt)
      expect_lt(max(abs(sp / a - 1)), 1e-8)
    }
  }
})

test_that("ptweedie() near p = 2 with large xi but moderate dispersion (p = 2.001, phi = 0.01)", {
  # xi = 1e5 here, but the dispersion is not small, so the saddlepoint
  # approximation is not accurate (relative error about 1e-6 to 1e-4); and the
  # tail between 1e-10 and 1e-5 used to keep the inversion value, wrong by 0.5%.
  # Reference values: the inversion integral evaluated in 240-bit arithmetic.
  p <- 2.001; phi <- 0.01
  expect_equal(ptweedie(1, mu = 1, phi = phi, power = p), 0.5133054316941893, tolerance = 1e-10)
  expect_equal(ptweedie(0.49999999999997202, mu = 1, phi = phi, power = p),
               3.1849241737973576e-10, tolerance = 1e-9)
  expect_equal(ptweedie(1.50000000000003, mu = 1, phi = phi, power = p, lower.tail = FALSE),
               5.9317066503315411e-06, tolerance = 1e-9)
})

test_that("ptweedie_saddle() does not return log-probabilities above zero far in the tails", {
  # The Mills ratio used to be computed as exp(log Phi(w) - log dnorm(w)),
  # which for |w| near 1e10 is a difference of two numbers near -5e19, so it
  # returned Inf, and ptweedie() then gave F = 1 for y far below the mean.
  expect_lt(ptweedie(1e-3, mu = 0.01, phi = 0.001, power = 8, log.p = TRUE), -1e18)
  expect_equal(ptweedie(1e-3, mu = 0.01, phi = 0.001, power = 8), 0)
  expect_equal(ptweedie(1e-5, mu = 0.01, phi = 0.1, power = 6), 0)
  # and the continued fraction agrees with the direct formula where both are accurate
  x <- c(3, 4, 6, 10, 20)
  expect_equal(tweedie:::mills_ratio(x),
               exp(pnorm(-x, log.p = TRUE) - dnorm(x, log = TRUE)), tolerance = 1e-13)
})
