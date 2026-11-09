# qtweedie(): the root-finding step (since 3.2.0)
#
# Before, when the starting value happened to be the quantile exactly (that
# is, ptweedie(start) == p), the starting value was stored in the wrong place
# (ans2[i] rather than ans[i]), and uniroot() was then called with an
# undefined upper bound (start.2), giving an error.

test_that("qtweedie() returns the starting value when it is the quantile exactly", {
  # For xi > 2, the starting value is the gamma quantile:
  #   qgamma(p, shape = 1/phi, rate = 1/(phi * mu)),
  # which is log(2) for p = 0.5, mu = 1, phi = 1.
  # Hitting  ptweedie(start) == p  exactly is very rare with the real
  # ptweedie(), so ptweedie() is replaced by a continuous, increasing
  # distribution function that equals 0.5 exactly on [0.6, 0.8] (which
  # contains log(2) = 0.693...), and increases elsewhere.
  start <- stats::qgamma(0.5, shape = 1, rate = 1)
  expect_equal(start, log(2))

  mock_ptweedie <- function(q, xi = NULL, mu, phi, power = NULL, ...) {
    ifelse(q < 0.6, 0.5 * q / 0.6,
           ifelse(q <= 0.8, 0.5,
                  0.5 + 0.5 * (1 - exp(-(q - 0.8)))))
  }
  local_mocked_bindings(ptweedie = mock_ptweedie, .package = "tweedie")

  # The single value
  expect_identical(qtweedie(0.5, xi = 2.5, mu = 1, phi = 1), start)

  # Within a vector, the value goes in the right place: here, p = 0 and p = 1
  # come first, so the position among the values solved for (i) differs from
  # the position in the output
  q <- qtweedie(c(0, 1, 0.5, 0.25, 0.75), xi = 2.5, mu = 1, phi = 1)
  expect_identical(q[1:3], c(0, Inf, start))
  # (the other values are the roots of the replacement distribution function)
  expect_equal(mock_ptweedie(q[4:5]), c(0.25, 0.75), tolerance = 1e-10)

  # And with names
  expect_identical(qtweedie(c(a = 0, b = 0.5), xi = 2.5, mu = 1, phi = 1),
                   c(a = 0, b = start))
})


test_that("qtweedie() inverts ptweedie()", {
  # Probabilities above Pr(Y = 0) when 1 < xi < 2 (below it, the quantile is 0)
  probs <- c(0.6, 0.7, 0.9, 0.99, 0.9999)
  for (xi in c(1.2, 1.5, 1.9, 2.5, 3, 4, 7)) {
    for (mu in c(0.3, 2)) {
      for (phi in c(0.2, 1, 5)) {
        if (xi < 2) {
          p0 <- dtweedie(0, xi = xi, mu = mu, phi = phi)   # Pr(Y = 0)
          pp <- probs[probs > p0 + 0.01]
          if (length(pp) == 0) next
        } else {
          pp <- c(0.001, 0.05, 0.3, 0.5, probs)
        }
        q <- suppressWarnings(qtweedie(pp, xi = xi, mu = mu, phi = phi))
        expect_true(all(is.finite(q)))
        expect_true(all(diff(q) > 0))          # (increasing in p)
        expect_equal(ptweedie(q, xi = xi, mu = mu, phi = phi), pp,
                     tolerance = 1e-8,
                     info = paste("xi =", xi, " mu =", mu, " phi =", phi))
      }
    }
  }
})


test_that("qtweedie() gives 0 below Pr(Y = 0) when 1 < xi < 2", {
  p0 <- dtweedie(0, xi = 1.5, mu = 1, phi = 1)        # exp(-2) = 0.135...
  q  <- qtweedie(c(p0 / 2, p0 * 0.999, 0.5), xi = 1.5, mu = 1, phi = 1)
  expect_identical(q[1:2], c(0, 0))
  expect_true(q[3] > 0)
})


test_that("qtweedie() with vector mu and phi matches element-by-element calls", {
  pp  <- c(0.2, 0.5, 0.8, 0.95)
  mu  <- c(0.5, 1, 2, 4)
  phi <- c(2, 1, 0.5, 0.25)
  for (xi in c(1.5, 3)) {
    q_vec <- qtweedie(pp, xi = xi, mu = mu, phi = phi)
    q_one <- mapply(function(p, m, f) qtweedie(p, xi = xi, mu = m, phi = f),
                    pp, mu, phi)
    expect_identical(q_vec, q_one)
  }
})
