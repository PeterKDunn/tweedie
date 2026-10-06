# Small tail probabilities for p > 2 are computed by integrating the density,
# since the Fourier inversion of the CDF has only ABSOLUTE accuracy (~1e-15).

test_that("tail integration matches the exact inverse Gaussian (p = 3, density by inversion)", {
  ti <- tweedie:::ptweedie_tail_integrate
  cases <- list(c(150, 1, 0), c(400, 1, 0), c(30, 0.1, 0), c(0.02, 0.1, 1), c(0.05, 1, 1))
  for (z in cases) {
    lt <- z[3] == 1
    r  <- ti(z[1], mu = 1, phi = z[2], power = 3, lower.tail = lt, IGexact = FALSE)
    ex <- statmod::pinvgauss(z[1], mean = 1, dispersion = z[2], lower.tail = lt, log.p = TRUE)
    expect_true(r$ok)
    expect_lt(abs(exp(r$logp - ex) - 1), 1e-10)   # relative error; tails down to ~1e-106
  }
})

test_that("ptweedie (p > 2) gives small tail probabilities with relative accuracy", {
  y <- c(20, 50, 100, 200)
  S <- ptweedie(y, mu = 1, phi = 1, power = 2.5, lower.tail = FALSE)
  expect_true(all(S > 0))                 # the inversion gave 0 or negative values for y >= 50
  expect_true(all(diff(log(S)) < 0))      # decreasing
  expect_equal(log(S), ptweedie(y, mu = 1, phi = 1, power = 2.5, lower.tail = FALSE, log.p = TRUE))
  # Where the inversion is still accurate (tail 1e-3 .. 1e-9), the two agree closely
  inv <- suppressWarnings(ptweedie_inversion(15, mu = 1, phi = 0.3, power = 6, lower.tail = FALSE))
  int <- exp(tweedie:::ptweedie_tail_integrate(15, mu = 1, phi = 0.3, power = 6, lower.tail = FALSE)$logp)
  expect_equal(int, inv, tolerance = 1e-7)
})

test_that("ptweedie is unchanged where the tail is not small", {
  y <- c(0.5, 1, 2, 5)
  expect_identical(ptweedie(y, mu = 1, phi = 1, power = 2.5),
                   suppressWarnings(ptweedie_inversion(y, mu = 1, phi = 1, power = 2.5)))
})
