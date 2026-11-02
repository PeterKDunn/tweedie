# The cap on the number of series terms (argument max.terms of ptweedie_series)

test_that("ptweedie_series() honours max.terms", {
  # lambda = 1e6 needs about 15 000 terms at y = mu
  phi <- 2 / 1e6
  full <- ptweedie_series(1, mu = 1, phi = phi, power = 1.5, details = TRUE)
  expect_gt(full$iterations, 1000)
  # a smaller cap truncates the sum, with a warning
  expect_warning(cut <- ptweedie_series(1, mu = 1, phi = phi, power = 1.5, max.terms = 1000),
                 "truncated")
  expect_lt(cut, full$cdf)
  # the default is unchanged, and ptweedie() works for 1 < p < 2
  expect_silent(ptweedie_series(1, mu = 1, phi = phi, power = 1.5))
  expect_equal(ptweedie(1, mu = 1, phi = phi, power = 1.5), full$cdf)
})
