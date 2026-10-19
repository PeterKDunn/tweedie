# The summation limits of ptweedie_series() must follow the terms
# Pr(N) * F_G (or Q_G), not just the Poisson weights; otherwise terms that
# matter are omitted far in either tail. The terms are log-concave in N, so
# the code finds their peak and sums outwards from it. References: 50-60 digit sums of the
# Poisson-gamma series (mpmath).

test_that("series keeps relative accuracy far in the upper tail", {
  # p = 1.5, mu = 1, phi = 4: Poisson mean 0.5, gamma shape 1, scale 2
  y   <- c(100, 200, 300)
  ref <- c(3.4136489462303752e-20, 2.4362963524509583e-40, 8.2746323426280853e-61)
  v   <- ptweedie_series(y, mu = 1, phi = 4, power = 1.5, lower.tail = FALSE)
  expect_lt(max(abs(v - ref) / ref), 1e-12)
  lv  <- ptweedie_series(y, mu = 1, phi = 4, power = 1.5, lower.tail = FALSE, log.p = TRUE)
  expect_lt(max(abs(lv - log(ref))), 1e-11)
})

test_that("series keeps relative accuracy far in the lower tail", {
  # p = 1.5, mu = 1, phi = 0.01: lambda = 200, gamma scale 0.005
  y   <- c(0.3, 0.5, 0.8)
  ref <- c(1.0093541776661044e-19, 2.8003157903556269e-9, 0.018533585771706545)
  v   <- ptweedie_series(y, mu = 1, phi = 0.01, power = 1.5)
  expect_lt(max(abs(v - ref) / ref), 1e-12)
})

test_that("series warns when truncated at the term cap", {
  # lambda = 1e11: about 5.6 million terms matter, more than the cap of 1e6
  expect_warning(ptweedie_series(1, mu = 1, phi = 1e-9, power = 1.99), "truncated")
})

test_that("series handles large lambda without truncation (lambda = 1e6)", {
  # Reference: direct sum of all terms N = 1, ..., 3e6 (and the atom)
  expect_silent(v <- ptweedie_series(0.1, mu = 0.1, phi = 0.001, power = 1.999, log.p = TRUE))
  expect_equal(v, -0.6847664295, tolerance = 1e-9)
})

test_that("series finds the terms that matter beyond N = 1e6 (far upper tail)", {
  # lambda = 55700, but the largest term is at N = 1188959.
  # Reference: direct sum of all terms N = 1, ..., 5e6
  d <- ptweedie_series(3e4, mu = 1000, phi = 0.01, power = 1.1, lower.tail = FALSE,
                       log.p = TRUE, details = TRUE)
  expect_equal(d$cdf, -3201730.8260874483, tolerance = 1e-12)
  expect_true(d$N.lo <= 1188959 && d$N.hi >= 1188959)
  expect_lt(d$iterations, 1e4)                 # only the terms that matter
})

test_that("series sums only the terms that matter, separately for each value", {
  d <- ptweedie_series(c(0.05, 1, 3), mu = 1, phi = 0.01, power = 1.05, details = TRUE)
  expect_length(d$iterations, 3)
  expect_true(all(d$N.lo <= d$N.hi))
  expect_lt(max(d$iterations), 400)
  # probabilities never exceed one, even when the sum rounds above it
  expect_true(all(ptweedie_series(c(1.5, 3, 10) * 1000, mu = 1000, phi = 0.001,
                                  power = 1.001, log.p = TRUE) <= 0))
})

test_that("series follows rising terms below the Poisson range (far lower tail, large lambda)", {
  # lambda ~ 105: the Poisson range starts near N = 40, but the terms that
  # matter are at small N; the point mass exp(-lambda) ~ 2e-46 is far smaller.
  # References: 50-digit series and 50-digit inversion, agreeing to 28+ digits.
  expect_equal(ptweedie_series(0.05, mu = 1, phi = 0.01, power = 1.05),
               5.3770368198611647265e-38, tolerance = 1e-12)
  expect_equal(ptweedie_series(0.05, mu = 1, phi = 0.01, power = 1.2),
               2.145118266660732397e-42, tolerance = 1e-12)
})
