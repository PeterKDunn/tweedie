# ptweedie(log.p = TRUE) must never return a positive log-probability.
# Rounding in the inversion can give a probability fractionally above one
# (by about 1e-14), whose log was returned as it stood before 3.2.1.

test_that("log-probabilities from ptweedie() are never positive", {
  # p > 2 (inversion): upper body of the distribution
  for (p in c(2.5, 4, 6)) {
    lp <- suppressWarnings(ptweedie(c(5, 10, 50, 100), mu = 1, phi = 0.01,
                                    power = p, log.p = TRUE))
    expect_true(all(lp <= 0))
  }
  # 1 < p < 2 with lambda > 1e9 (inversion)
  lp <- suppressWarnings(ptweedie(c(1.2, 2, 5), mu = 1, phi = 1e-10,
                                  power = 1.5, log.p = TRUE))
  expect_true(all(lp <= 0))
  # and the linear-scale values agree: exp(log F) = F
  F  <- suppressWarnings(ptweedie(c(5, 10), mu = 1, phi = 0.01, power = 2.5))
  lp <- suppressWarnings(ptweedie(c(5, 10), mu = 1, phi = 0.01, power = 2.5, log.p = TRUE))
  expect_equal(exp(lp), F, tolerance = 1e-13)
})
