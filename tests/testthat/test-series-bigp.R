# For p > 2 the series for the density alternates in sign; when its terms
# peak at a large index it cannot be summed accurately in double precision.
# It used to return 0, Inf, absurd values (e.g. 1e239), or hang; it now
# returns NaN (with a warning), and dtweedie() falls back to other methods.

test_that("dtweedie_series returns NaN with a warning when it cannot be summed", {
  expect_warning(f <- dtweedie_series(0.005, power = 2.05, mu = 0.01, phi = 0.01),
                 "cannot be summed accurately")
  expect_true(is.nan(f))
  # Previously returned 7.98e+102 (true density ~ 3.99)
  expect_warning(f <- dtweedie_series(1, power = 2.5, mu = 1, phi = 0.01),
                 "cannot be summed accurately")
  expect_true(is.nan(f))
})

test_that("dtweedie_series does not hang for a huge peak index", {
  # Peak index ~ 7.7e24: the search for summation limits used to loop forever
  expect_warning(f <- dtweedie_series(0.01, power = 15, mu = 1, phi = 10),
                 "cannot be summed accurately")
  expect_true(is.nan(f))
})

test_that("dtweedie_series is unchanged where it can be summed", {
  f <- as.numeric(dtweedie_series(c(1, 2, 5), power = 3.5, mu = 1, phi = 1))
  expect_true(all(is.finite(f)))
  expect_equal(f, dtweedie_inversion(c(1, 2, 5), mu = 1, phi = 1, power = 3.5), tolerance = 1e-8)
})

test_that("dtweedie for p > 10 no longer hangs or returns Inf", {
  # Near the mean the density is ~ 1 / (sd * sqrt(2 pi)), sd = sqrt(phi mu^p)
  f <- dtweedie(0.05, mu = 0.05, phi = 0.05, power = 12)
  expect_equal(f, dtweedie_inversion(0.05, mu = 0.05, phi = 0.05, power = 12), tolerance = 1e-8)
  expect_equal(f, 1 / (sqrt(0.05 * 0.05^12) * sqrt(2 * pi)), tolerance = 1e-2)
  expect_equal(dtweedie(0.1, mu = 0.1, phi = 10, power = 11),
               dtweedie_inversion(0.1, mu = 0.1, phi = 10, power = 11), tolerance = 1e-8)
  # Extremely far in the lower tail: density underflows to 0 (was Inf / hang)
  expect_equal(dtweedie(0.001, mu = 1, phi = 0.05, power = 12), 0)
  expect_equal(dtweedie(0.01,  mu = 1, phi = 10,   power = 15), 0)
})
