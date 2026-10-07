# (as.numeric() because dtweedie_series() returns an array with a dim attribute)
# For p = 1, Y = phi * N with N ~ Poisson(mu / phi), so Y lives on the lattice
# 0, phi, 2 phi, ...; the "density" is the probability Pr(Y = y) = Pr(N = y/phi).
# dtweedie() and dtweedie_series() must agree, and match ptweedie().

test_that("p = 1: dtweedie_series() agrees with dtweedie() when phi != 1", {
  for (phi in c(0.1, 0.5, 2)) {
    y <- phi * (0:10)
    expect_equal(as.numeric(dtweedie_series(y, mu = 1, phi = phi, power = 1)),
                 dtweedie(y, mu = 1, phi = phi, power = 1))
  }
})

test_that("p = 1: the density is the Poisson probability Pr(N = y/phi)", {
  y <- c(0, 0.1, 0.3, 0.7)           # note 0.3/0.1 is not exactly 3 in floating point
  expect_equal(as.numeric(dtweedie_series(y, mu = 1, phi = 0.1, power = 1)),
               dpois(c(0, 1, 3, 7), lambda = 10))
  expect_equal(dtweedie(y, mu = 1, phi = 0.1, power = 1),
               dpois(c(0, 1, 3, 7), lambda = 10))
})

test_that("p = 1: the densities sum to one, and match the jumps in ptweedie()", {
  phi <- 0.5
  y <- phi * (0:60)
  d <- dtweedie(y, mu = 3, phi = phi, power = 1)
  expect_equal(sum(d), 1, tolerance = 1e-12)
  expect_equal(cumsum(d), ptweedie(y, mu = 3, phi = phi, power = 1), tolerance = 1e-12)
})

test_that("p = 1, phi = 1: the ordinary Poisson distribution", {
  expect_equal(as.numeric(dtweedie_series(0:8, mu = 2.5, phi = 1, power = 1)), dpois(0:8, 2.5))
})
