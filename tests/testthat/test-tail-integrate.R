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

test_that("upper tails for p > 2 (Gauss-Laguerre path) match the exact inverse Gaussian", {
  # The general algorithm (IGexact = FALSE) against statmod::pinvgauss, for
  # upper-tail probabilities from 1e-10 down to about 1e-290
  cases <- expand.grid(q = c(5, 10, 20, 50, 100, 150, 300), mu = c(0.5, 1.4, 5), phi = c(0.1, 0.74, 2))
  cases$exact <- statmod::pinvgauss(cases$q, cases$mu, dispersion = cases$phi, lower.tail = FALSE)
  cases <- cases[cases$exact < 1e-10 & cases$exact > 1e-300, ]
  expect_gt(nrow(cases), 20)
  for (i in seq_len(nrow(cases))) {
    r <- tweedie:::ptweedie_tail_integrate(cases$q[i], cases$mu[i], cases$phi[i], 3,
                                           lower.tail = FALSE, IGexact = FALSE)
    expect_true(r$ok)
    expect_lt(abs(expm1(r$logp - log(cases$exact[i]))), 1e-10)
  }
})

test_that("upper tails for p > 2 need few density evaluations", {
  # Gauss-Laguerre with 10 and 20 nodes: 30 evaluations, plus one at q
  n <- 0
  trace("dtweedie_inversion", function() n <<- n + length(get("y", parent.frame())),
        print = FALSE, where = asNamespace("tweedie"))
  on.exit(untrace("dtweedie_inversion", where = asNamespace("tweedie")))
  r <- tweedie:::ptweedie_tail_integrate(150, 1.4, 0.74, 3, lower.tail = FALSE, IGexact = FALSE)
  expect_true(r$ok)
  expect_lte(n, 31)
})

test_that("lower tails for p > 2 (Gauss-Laguerre path) match the exact inverse Gaussian", {
  # The general algorithm (IGexact = FALSE) against statmod::pinvgauss, for
  # lower-tail probabilities from 1e-5 down to about 1e-290
  cases <- expand.grid(q = c(0.002, 0.005, 0.01, 0.03, 0.1, 0.3), mu = c(0.5, 1.4, 5), phi = c(0.1, 0.74, 2))
  cases$exact <- statmod::pinvgauss(cases$q, cases$mu, dispersion = cases$phi)
  cases <- cases[cases$exact < 1e-5 & cases$exact > 1e-290, ]
  expect_gt(nrow(cases), 15)
  for (i in seq_len(nrow(cases))) {
    r <- tweedie:::ptweedie_tail_integrate(cases$q[i], cases$mu[i], cases$phi[i], 3,
                                           lower.tail = TRUE, IGexact = FALSE)
    expect_true(r$ok)
    expect_lt(abs(expm1(r$logp - log(cases$exact[i]))), 1e-10)
  }
})

test_that("lower tails for p > 2 need few density evaluations", {
  n <- 0
  trace("dtweedie_inversion", function() n <<- n + length(get("y", parent.frame())),
        print = FALSE, where = asNamespace("tweedie"))
  on.exit(untrace("dtweedie_inversion", where = asNamespace("tweedie")))
  r <- tweedie:::ptweedie_tail_integrate(0.02, 1, 0.5, 3, lower.tail = TRUE, IGexact = FALSE)
  expect_true(r$ok)
  expect_lte(n, 31)
})

test_that("ptweedie (p > 2) has relative accuracy for tails between 1e-10 and 1e-5", {
  # References: the piecewise adaptive integration of the density (an
  # independent route). The inversion alone is in error by 1e-7 and 5e-8 here.
  S <- ptweedie(25, mu = 1, phi = 1, power = 2.5, lower.tail = FALSE)
  expect_lt(abs(S / 2.8982324446682487e-09 - 1), 1e-12)
  F <- ptweedie(0.005, mu = 1, phi = 1, power = 2.5)
  expect_lt(abs(F / 7.3715010314434954e-09 - 1), 1e-12)
  # and the log scale agrees
  expect_equal(ptweedie(25, mu = 1, phi = 1, power = 2.5, lower.tail = FALSE, log.p = TRUE), log(S))
})
