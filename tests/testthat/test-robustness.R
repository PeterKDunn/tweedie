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
  expect_equal(dtweedie(1e-10, mu = 1, phi = 1, power = 50), 0)
  expect_equal(dtweedie(1e-10, mu = 1e3, phi = 1e3, power = 50), 0)
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
