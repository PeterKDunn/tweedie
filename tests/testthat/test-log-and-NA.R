# dtweedie(log = TRUE), and NA in ptweedie() (since 3.2.0)

test_that("dtweedie(log = TRUE) agrees with log(dtweedie()) where the density does not underflow", {
  y <- c(0, 0.01, 0.5, 1, 2, 7)
  for (p in c(1, 1.05, 1.5, 1.95, 2, 2.05, 2.5, 3, 4, 8)) {
    for (phi in c(0.05, 1)) {
      L <- suppressWarnings(dtweedie(y, mu = 1.3, phi = phi, power = p, log = TRUE))  # (p = 1: non-integer y warns)
      f <- suppressWarnings(dtweedie(y, mu = 1.3, phi = phi, power = p))
      ok <- f > 0
      expect_equal(L[ok], log(f[ok]), tolerance = 1e-13)
      # where the density underflows (or is zero), the log is below log(1e-308)
      expect_true(all(L[!ok] < -700))
    }
  }
})

test_that("dtweedie(log = TRUE) stays finite where the density underflows", {
  # References: the series summed in Rmpfr (on the log scale)
  expect_equal(dtweedie(50,  mu = 1, phi = 0.01, power = 1.5, log = TRUE), -7373.123378554057712, tolerance = 1e-14)
  expect_equal(dtweedie(500, mu = 1, phi = 0.01, power = 2.5, log = TRUE), -33145.69159941934986, tolerance = 1e-14)
  expect_equal(dtweedie(50,  mu = 1, phi = 0.01, power = 4,   log = TRUE), -1622.800418613131322, tolerance = 1e-14)
  expect_equal(dtweedie(50,  mu = 1, phi = 0.01, power = 1.01, log = TRUE), -14427.31681252414503, tolerance = 1e-14)
  expect_equal(dtweedie(50,  mu = 1, phi = 0.01, power = 1.5), 0)
  # y = 0 (the point mass, 1 < p < 2), NA, Inf, names
  expect_equal(dtweedie(c(a = 0, b = NA, c = Inf), mu = 1, phi = 1, power = 1.5, log = TRUE), 
               c(a = -2, b = NA, c = -Inf))
  expect_equal(dtweedie(0, mu = 1, phi = 1, power = 2.5, log = TRUE), -Inf)
})

test_that("the log-likelihood does not become -Inf for observations far in the tails", {
  # One observation far in the tail: the density underflows, but its log does not
  y   <- c(1, 1.1, 0.9, 50)
  ll  <- tweedie:::dtweedie_logl(phi = 0.01, y = y, mu = 1, power = 1.5)
  expect_true(is.finite(ll))
  expect_equal(ll, sum(dtweedie(y, mu = 1, phi = 0.01, power = 1.5, log = TRUE)))
})

test_that("ptweedie() gives NA for NA, and the other values as usual", {
  for (p in c(1.5, 2, 2.5, 3)) {
    expect_equal(ptweedie(c(1, NA, 2), mu = 1, phi = 1, power = p),
                 c(ptweedie(1, mu = 1, phi = 1, power = p), NA, ptweedie(2, mu = 1, phi = 1, power = p)))
  }
  expect_equal(ptweedie(c(a = 1, b = NA), mu = c(1, 2), phi = 1, power = 2.5),
               c(a = ptweedie(1, mu = 1, phi = 1, power = 2.5), b = NA))
  expect_equal(ptweedie(c(NA, NA), mu = 1, phi = 1, power = 1.5), c(NA_real_, NA_real_))
  expect_equal(ptweedie(c(1, NA), mu = 1, phi = 1, power = 1.5, lower.tail = FALSE, log.p = TRUE),
               c(ptweedie(1, mu = 1, phi = 1, power = 1.5, lower.tail = FALSE, log.p = TRUE), NA))
})
