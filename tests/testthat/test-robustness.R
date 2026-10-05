# Inputs that previously returned NaN, raised an error, or returned an
# impossible value from the user-level functions.

test_that("ptweedie upper tail at y <= 0 is 1 when there is no mass there", {
  for (p in c(2.5, 3.5, 5)) {
    expect_equal(ptweedie(0,  mu = 1, phi = 1, power = p, lower.tail = FALSE), 1)
    expect_equal(ptweedie(0,  mu = 1, phi = 1, power = p, lower.tail = TRUE),  0)
    expect_equal(ptweedie(0,  mu = 1, phi = 1, power = p, lower.tail = FALSE, log.p = TRUE), 0)
  }
  expect_equal(ptweedie(-1, mu = 1, phi = 1, power = 1.5, lower.tail = FALSE), 1)
  expect_equal(ptweedie(-1, mu = 1, phi = 1, power = 1.5), 0)
})

test_that("dtweedie with very large p and tiny y returns 0, not an error", {
  # y^(2-p) overflows, so the deviance is Inf and the density is 0
  expect_silent(f1 <- dtweedie(1e-10, mu = 1, phi = 1, power = 50))
  expect_silent(f2 <- dtweedie(1e-10, mu = 1e3, phi = 1e3, power = 50))
  expect_equal(c(f1, f2), c(0, 0))
})

test_that("ptweedie far in the upper tail (1 < p < 2) is not NaN", {
  # The inversion returned NaN here; the series (positive terms) is used instead
  expect_equal(ptweedie(1e3, mu = 1e-3, phi = 10, power = 1.01, lower.tail = FALSE), 0)
  expect_equal(ptweedie(1e6, mu = 1e-3, phi = 1e3, power = 1.01), 1)
  expect_false(is.nan(ptweedie(1e3, mu = 1e-3, phi = 1e3, power = 1.5, lower.tail = FALSE)))
})

test_that("ptweedie no longer returns 0.5 for tiny y when the inversion fails (1 < p < 2)", {
  # P(Y = 0) = exp(-1) here, so F(1e-10) is a little above 0.368
  F <- ptweedie(1e-10, mu = 1, phi = 10, power = 1.9)
  expect_equal(F, as.numeric(ptweedie_series(1e-10, mu = 1, phi = 10, power = 1.9)), tolerance = 1e-10)
  expect_gt(F, exp(-1))
  expect_lt(ptweedie(1e-10, mu = 1e3, phi = 1, power = 1.5), 1e-20)
})

test_that("ptweedie_series handles vector mu and phi", {
  v <- ptweedie_series(c(0.5, 1, 2), power = 1.5, mu = c(1, 2, 3), phi = c(1, 0.5, 2))
  s <- mapply(function(q, m, f) ptweedie_series(q, power = 1.5, mu = m, phi = f),
              c(0.5, 1, 2), c(1, 2, 3), c(1, 0.5, 2))
  expect_equal(v, s)
  expect_silent(ptweedie_series(c(0.5, 1, 2), power = 1.5, mu = c(1, 1, 1), phi = 1))
})

test_that("p near 1 with small phi: no longer converges falsely (was silently wrong)", {
  # The convergence test used the integrand's value at a zero (~0 by
  # construction) instead of its amplitude, and stopped after 4 regions with
  # F wrong by ~0.02 and exitstatus = 0. The series is the reference here.
  for (z in list(c(0.3, 0.001, 1.05), c(0.05, 0.01, 1.01), c(0.3, 0.01, 1.2))) {
    o <- suppressWarnings(ptweedie_inversion(z[1], mu = 1, phi = z[2], power = z[3], details = TRUE))
    s <- as.numeric(ptweedie_series(z[1], mu = 1, phi = z[2], power = z[3]))
    expect_true(abs(o$cdf - s) < 1e-6 || o$exitstatus == 1L)   # accurate, or else flagged
  }
  # and this one, from an earlier fix, must still converge quickly
  o <- ptweedie_inversion(1, mu = 5, phi = 2, power = 1.01, details = TRUE)
  expect_lt(o$regions, 50)
  expect_equal(o$cdf, as.numeric(ptweedie_series(1, mu = 5, phi = 2, power = 1.01)), tolerance = 1e-8)
})

test_that("non-convergence is reported by a warning, not only through exitstatus", {
  expect_warning(ptweedie(0.999, mu = 1, phi = 1e-6, power = 2.5),
                 "did not reach the target accuracy")
  expect_warning(o <- ptweedie_inversion(0.999, mu = 1, phi = 1e-6, power = 2.5, details = TRUE),
                 "did not reach the target accuracy")
  expect_equal(o$exitstatus, 1L)
  # no warning when everything converges
  expect_silent(ptweedie(c(0.5, 1, 2), mu = 1, phi = 1, power = 2.5))
  expect_silent(dtweedie(c(0.5, 1, 2), mu = 1, phi = 1, power = 1.5))
})
