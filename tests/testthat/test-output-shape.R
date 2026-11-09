# The length and form of the output of dtweedie(), ptweedie(), qtweedie() and
# rtweedie() match those of R's own functions for distributions (since 3.2.0):
#
# * d-, p- and q-functions: a double vector the length of the first input,
#   with all the attributes of the first input (dim, dimnames, names, tsp/class,
#   ...); NA in gives NA out; a zero-length input gives a zero-length output.
# * r-functions: a plain double vector of length n (no attributes); if
#   length(n) > 1, its length is used; non-integer n is truncated; n = 0 gives
#   a zero-length result.
#
# Each case is compared directly with the corresponding base-R gamma function,
# so the tests check consistency with R itself.

# Values of xi covering each regime: Poisson (phi = 1), Poisson-gamma, gamma,
# p > 2 (series/inversion), inverse Gaussian, and large p
xi_values <- c(1, 1.5, 2, 2.5, 3, 4)

# Inputs of different forms (all valid as y, q and p: values in (0, 1))
shape_inputs <- function() {
  list(
    vector        = c(0.2, 0.4, 0.6),
    named         = c(a = 0.2, b = 0.4, c = 0.6),
    matrix        = matrix(c(0.1, 0.3, 0.5, 0.7), nrow = 2,
                           dimnames = list(c("r1", "r2"), c("c1", "c2"))),
    array_1d      = array(c(0.2, 0.4, 0.6)),
    array_3d      = array(seq(0.05, 0.95, length = 8), dim = c(2, 2, 2)),
    ts            = stats::ts(c(0.2, 0.4, 0.6), start = 2000),
    other_attr    = structure(c(0.2, 0.4), myattr = "kept"),
    scalar        = 0.5
  )
}

# Same attributes and type as the base-R function gives (values are not compared)
expect_same_form <- function(object, reference) {
  expect_identical(typeof(object), "double")
  expect_identical(length(object), length(reference))
  expect_identical(attributes(object), attributes(reference))
}

dfun <- function(x, xi) {
  # (xi = 1: non-integer y gives a warning, as for dpois())
  suppressWarnings(dtweedie(x, xi = xi, mu = 1, phi = 1))
}


test_that("dtweedie() output has the form of the first input, as for dgamma()", {
  for (xi in xi_values) {
    for (x in shape_inputs()) {
      expect_same_form(dfun(x, xi), dgamma(x, shape = 2))
      expect_same_form(suppressWarnings(dtweedie(x, xi = xi, mu = 1, phi = 1, log = TRUE)),
                       dgamma(x, shape = 2, log = TRUE))
    }
  }
})


test_that("ptweedie() output has the form of the first input, as for pgamma()", {
  for (xi in xi_values) {
    for (x in shape_inputs()) {
      expect_same_form(ptweedie(x, xi = xi, mu = 1, phi = 1),
                       pgamma(x, shape = 2))
      expect_same_form(ptweedie(x, xi = xi, mu = 1, phi = 1, lower.tail = FALSE, log.p = TRUE),
                       pgamma(x, shape = 2, lower.tail = FALSE, log.p = TRUE))
    }
  }
})


test_that("qtweedie() output has the form of the first input, as for qgamma()", {
  for (xi in xi_values) {
    for (x in shape_inputs()) {
      expect_same_form(suppressWarnings(qtweedie(x, xi = xi, mu = 1, phi = 1)),
                       qgamma(x, shape = 2))
    }
  }
})


test_that("the values themselves are unchanged by the form of the input", {
  m <- matrix(c(0.1, 0.3, 0.5, 0.7), nrow = 2)
  for (xi in c(1.5, 2.5)) {
    expect_equal(as.vector(dtweedie(m, xi = xi, mu = 1, phi = 1)),
                 dtweedie(as.vector(m), xi = xi, mu = 1, phi = 1))
    expect_equal(as.vector(ptweedie(m, xi = xi, mu = 1, phi = 1)),
                 ptweedie(as.vector(m), xi = xi, mu = 1, phi = 1))
    expect_equal(as.vector(qtweedie(m, xi = xi, mu = 1, phi = 1)),
                 qtweedie(as.vector(m), xi = xi, mu = 1, phi = 1))
  }
  # A vector mu (same length as the first input) with a matrix first input
  mu <- c(0.5, 1, 2, 4)
  expect_equal(as.vector(dtweedie(m, xi = 1.5, mu = mu, phi = 1)),
               dtweedie(as.vector(m), xi = 1.5, mu = mu, phi = 1))
  expect_identical(dim(dtweedie(m, xi = 1.5, mu = mu, phi = 1)), c(2L, 2L))
})


test_that("integer and logical first inputs give double results, as for dgamma()", {
  for (xi in c(1.5, 3)) {
    expect_same_form(dtweedie(1:3, xi = xi, mu = 1, phi = 1), dgamma(1:3, shape = 2))
    expect_same_form(ptweedie(1:3, xi = xi, mu = 1, phi = 1), pgamma(1:3, shape = 2))
    expect_same_form(dtweedie(c(TRUE, FALSE), xi = xi, mu = 1, phi = 1),
                     dgamma(c(TRUE, FALSE), shape = 2))
  }
})


test_that("zero-length inputs give zero-length results", {
  for (xi in xi_values) {
    expect_identical(dtweedie(numeric(0), xi = xi, mu = 1, phi = 1), numeric(0))
    expect_identical(ptweedie(numeric(0), xi = xi, mu = 1, phi = 1), numeric(0))
    expect_identical(qtweedie(numeric(0), xi = xi, mu = 1, phi = 1), numeric(0))
  }
})


test_that("NA in gives NA out, keeping the form of the first input", {
  x <- c(a = 0.2, b = NA, c = 0.6)
  m <- matrix(c(0.1, NA, 0.5, 0.7), nrow = 2)
  for (xi in xi_values) {
    d <- dfun(x, xi)
    expect_same_form(d, dgamma(x, shape = 2))
    expect_true(is.na(d[["b"]]))
    expect_false(anyNA(d[c("a", "c")]))

    pr <- ptweedie(x, xi = xi, mu = 1, phi = 1)
    expect_same_form(pr, pgamma(x, shape = 2))
    expect_true(is.na(pr[["b"]]))
    expect_equal(pr[c("a", "c")], ptweedie(x[c("a", "c")], xi = xi, mu = 1, phi = 1))
    expect_same_form(ptweedie(m, xi = xi, mu = 1, phi = 1), pgamma(m, shape = 2))

    qu <- suppressWarnings(qtweedie(x, xi = xi, mu = 1, phi = 1))
    expect_same_form(qu, qgamma(x, shape = 2))
    expect_true(is.na(qu[["b"]]))
    expect_equal(qu[c("a", "c")],
                 suppressWarnings(qtweedie(x[c("a", "c")], xi = xi, mu = 1, phi = 1)))
    expect_same_form(suppressWarnings(qtweedie(m, xi = xi, mu = 1, phi = 1)),
                     qgamma(m, shape = 2))
  }
  # All values missing
  expect_identical(qtweedie(c(NA_real_, NA_real_), xi = 1.5, mu = 1, phi = 1),
                   c(NA_real_, NA_real_))
  # NA with a vector mu: the remaining mu values stay matched to their p values
  mu <- c(1, 2, 3)
  expect_equal(qtweedie(c(0.3, NA, 0.7), xi = 1.5, mu = mu, phi = 1)[c(1, 3)],
               qtweedie(c(0.3, 0.7), xi = 1.5, mu = mu[c(1, 3)], phi = 1))
})


test_that("qtweedie() handles probabilities of exactly 0 and 1", {
  for (xi in c(1.5, 2.5)) {
    expect_identical(qtweedie(c(0, 1), xi = xi, mu = 1, phi = 1), c(0, Inf))
    expect_identical(qtweedie(c(a = 0, b = 1), xi = xi, mu = 1, phi = 1), c(a = 0, b = Inf))
    q3 <- qtweedie(c(0, 0.5, 1), xi = xi, mu = 1, phi = 1)
    expect_identical(q3[c(1, 3)], c(0, Inf))
    expect_true(q3[2] > 0 && is.finite(q3[2]))
  }
})


test_that("rtweedie() gives a plain double vector of length n, as for rgamma()", {
  set.seed(1)
  for (xi in xi_values) {
    for (n in c(1, 3, 10)) {
      r <- rtweedie(n, xi = xi, mu = 1, phi = 1)
      expect_same_form(r, rgamma(n, shape = 2))
      expect_null(attributes(r))
      expect_null(dim(r))
    }
  }
})


test_that("rtweedie() interprets n as rgamma() does", {
  set.seed(1)
  for (xi in xi_values) {
    # length(n) > 1: the length is the number required
    expect_length(rtweedie(c(5, 6, 7), xi = xi, mu = 1, phi = 1), 3)
    expect_length(rtweedie(c(a = 9), xi = xi, mu = 1, phi = 1), 9)
    # non-integer n is truncated
    expect_length(rtweedie(2.7, xi = xi, mu = 1, phi = 1), length(rgamma(2.7, shape = 2)))
    # n = 0, or a zero-length n, gives a zero-length result
    expect_identical(rtweedie(0, xi = xi, mu = 1, phi = 1), numeric(0))
    expect_identical(rtweedie(numeric(0), xi = xi, mu = 1, phi = 1), numeric(0))
    # mu and phi longer than n: the first n are used
    expect_length(rtweedie(2, xi = xi, mu = 1:5, phi = 1), 2)
    # invalid n
    expect_error(rtweedie(-1, xi = xi, mu = 1, phi = 1))
    expect_error(rtweedie(NA, xi = xi, mu = 1, phi = 1))
    expect_error(rtweedie(Inf, xi = xi, mu = 1, phi = 1))
  }
})


test_that("rtweedie() values are valid for each regime", {
  set.seed(1)
  for (xi in xi_values) {
    r <- rtweedie(50, xi = xi, mu = 2, phi = 1)
    expect_false(anyNA(r))
    expect_true(all(r >= 0))
    if (xi >= 2) expect_true(all(r > 0))
  }
  # xi = 1, phi = 2: values are multiples of phi
  r <- rtweedie(20, xi = 1, mu = 2, phi = 2)
  expect_equal(r %% 2, rep(0, 20))
})


test_that("tweedie_profile() works with rtweedie() data and with a 1-d array response", {
  set.seed(96)
  y <- rtweedie(100, xi = 1.5, mu = 1, phi = 1)
  # (capture.output() hides the progress bar)
  utils::capture.output(
    out <- tweedie_profile(y ~ 1, xi.vec = seq(1.3, 1.7, by = 0.1),
                           do.plot = FALSE, verbose = 0) )
  expect_true(is.finite(out$xi.max))
  expect_true(out$xi.max > 1 && out$xi.max < 2)

  # A 1-d array response (as, e.g., from tapply()) gives the same answer as a vector
  y_1d <- array(y)
  utils::capture.output(
    out_1d <- tweedie_profile(y_1d ~ 1, xi.vec = seq(1.3, 1.7, by = 0.1),
                              do.plot = FALSE, verbose = 0) )
  expect_equal(out_1d$xi.max, out$xi.max)
  expect_equal(out_1d$L, out$L)
})
