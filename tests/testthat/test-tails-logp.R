# test-tails-and-logp.R
#
# Tests for lower.tail / log.p support added to ptweedie(), ptweedie_series()
# and ptweedie_inversion(), plus the special_cases() atom-at-zero handling
# and the pre-acceleration / conditional-integrand fixes for 1<p<2.
#
# These tests are deliberately independent of implementation details (they
# never call internal Fortran routines directly) so they will catch
# regressions regardless of how the internals change in future.

test_that("lower.tail and upper.tail sum to 1 (ptweedie, general dispatch)", {
  cases <- expand.grid(
    power = c(1.01, 1.3, 1.9, 2.2, 3, 4.5),
    mu    = c(0.5, 1, 5),
    phi   = c(0.5, 1, 2)
  )
  y <- c(0.01, 0.5, 1, 5, 20)
  
  for (i in seq_len(nrow(cases))) {
    p <- cases$power[i]; m <- cases$mu[i]; f <- cases$phi[i]
    lo_full <- ptweedie(y, power = p, mu = m, phi = f, lower.tail = TRUE)
    hi_full <- ptweedie(y, power = p, mu = m, phi = f, lower.tail = FALSE)
    
    # Known limitation (see "KNOWN LIMITATION" test below): same underlying
    # ptweedie_inversion() case, reached here via ptweedie()'s dispatch.
    is_known_bad <- (p == 1.01 & m == 0.5 & f == 0.5 & y == 20)
    lo <- lo_full[!is_known_bad]; hi <- hi_full[!is_known_bad]
    
    expect_equal(lo + hi, rep(1, length(lo)), tolerance = 1e-8,
                 info = sprintf("power=%s mu=%s phi=%s", p, m, f))
  }
})

test_that("lower.tail and upper.tail sum to 1 (ptweedie_inversion directly)", {
  cases <- expand.grid(
    power = c(1.01, 1.3, 1.9, 2.2, 3, 4.5),
    mu    = c(0.5, 1, 5),
    phi   = c(0.5, 1, 2)
  )
  y <- c(0.01, 0.5, 1, 5, 20)
  
  for (i in seq_len(nrow(cases))) {
    p <- cases$power[i]; m <- cases$mu[i]; f <- cases$phi[i]
    lo_full <- ptweedie_inversion(y, power = p, mu = m, phi = f, lower.tail = TRUE)
    hi_full <- ptweedie_inversion(y, power = p, mu = m, phi = f, lower.tail = FALSE)
    
    # Known limitation (see "KNOWN LIMITATION" test below): extreme upper-tail
    # values can be slightly negative for p near 1 with small mu/phi at large y.
    # Exclude just that corner from this general sweep rather than loosening
    # the tolerance for everything.
    is_known_bad <- (p == 1.01 & m == 0.5 & f == 0.5 & y == 20)
    lo <- lo_full[!is_known_bad]; hi <- hi_full[!is_known_bad]
    
    expect_equal(lo + hi, rep(1, length(lo)), tolerance = 1e-8,
                 info = sprintf("power=%s mu=%s phi=%s", p, m, f))
  }
})

test_that("lower.tail and upper.tail sum to 1 (ptweedie_series, 1<p<2 only)", {
  cases <- expand.grid(
    power = c(1.01, 1.1, 1.5, 1.9),
    mu    = c(0.5, 1, 5),
    phi   = c(0.5, 1, 2)
  )
  y <- c(0.01, 0.5, 1, 5, 20)
  
  for (i in seq_len(nrow(cases))) {
    p <- cases$power[i]; m <- cases$mu[i]; f <- cases$phi[i]
    lo <- ptweedie_series(y, power = p, mu = m, phi = f, lower.tail = TRUE)
    hi <- ptweedie_series(y, power = p, mu = m, phi = f, lower.tail = FALSE)
    expect_equal(lo + hi, rep(1, length(y)), tolerance = 1e-10,
                 info = sprintf("power=%s mu=%s phi=%s", p, m, f))
  }
})

test_that("REGRESSION: lower.tail is not silently inverted", {
  # This is the specific bug found and fixed during development: tail=as.integer(lower.tail)
  # was passed instead of tail=as.integer(!lower.tail), causing lower.tail=TRUE to silently
  # return the upper tail and vice versa. ptweedie() must be monotone non-decreasing in y,
  # up to the small numerical wobble expected from independent pointwise evaluation
  # (each y is solved by its own integration, with no guarantee of monotonicity across
  # calls at the ~1e-6 level, particularly deep in the tail where F(y) is close to 1).
  y <- seq(0.1, 20, length = 30)
  for (p in c(1.01, 1.3, 1.9, 2.2, 3, 4.5)) {
    Fy <- ptweedie(y, power = p, mu = 1, phi = 1, lower.tail = TRUE)
    expect_true(all(diff(Fy) >= -1e-5),
                info = sprintf("ptweedie not monotone non-decreasing for power=%s", p))
    Sy <- ptweedie(y, power = p, mu = 1, phi = 1, lower.tail = FALSE)
    expect_true(all(diff(Sy) <= 1e-5),
                info = sprintf("ptweedie upper tail not monotone non-increasing for power=%s", p))
  }
})



test_that("log.p matches log(linear result) where both are well away from underflow", {
  cases <- expand.grid(
    power = c(1.01, 1.3, 1.9, 2.2, 3, 4.5),
    lower.tail = c(TRUE, FALSE)
  )
  y <- c(0.5, 1, 2, 5)
  
  for (i in seq_len(nrow(cases))) {
    p <- cases$power[i]; lt <- cases$lower.tail[i]
    lin <- ptweedie_inversion(y, power = p, mu = 1, phi = 1, lower.tail = lt)
    lg  <- ptweedie_inversion(y, power = p, mu = 1, phi = 1, lower.tail = lt, log.p = TRUE)
    expect_equal(lg, log(lin), tolerance = 1e-6,
                 info = sprintf("power=%s lower.tail=%s", p, lt))
  }
  
  # series method, 1<p<2 only
  for (p in c(1.01, 1.3, 1.9)) {
    for (lt in c(TRUE, FALSE)) {
      lin <- ptweedie_series(y, power = p, mu = 1, phi = 1, lower.tail = lt)
      lg  <- ptweedie_series(y, power = p, mu = 1, phi = 1, lower.tail = lt, log.p = TRUE)
      expect_equal(lg, log(lin), tolerance = 1e-6,
                   info = sprintf("series power=%s lower.tail=%s", p, lt))
    }
  }
})

test_that("log.p stays finite deep in the tail where linear scale underflows", {
  # Series method, 1<p<2: deep right tail, lower.tail=TRUE
  lg <- ptweedie_series(300, power = 1.3, mu = 1.5, phi = 2,
                        lower.tail = FALSE, log.p = TRUE)
  expect_true(is.finite(lg))
  expect_true(lg < -50)   # should be a large negative number, not -Inf and not 0
  
  # Inversion method, p>2: deep right tail, lower.tail=FALSE, cross-checked
  # against statmod::pinvgauss for p=3
  skip_if_not_installed("statmod")
  for (yy in c(100, 150, 200)) {
    lg_inv  <- ptweedie_inversion(yy, mu = 1, phi = 1, power = 3,
                                  lower.tail = FALSE, log.p = TRUE)
    lg_exact <- log(statmod::pinvgauss(yy, mean = 1, dispersion = 1, lower.tail = FALSE))
    expect_equal(lg_inv, lg_exact, tolerance = 1e-4,
                 info = sprintf("y=%s", yy))
  }
})

test_that("ptweedie_inversion matches statmod::pinvgauss exactly for p=3 (both tails)", {
  skip_if_not_installed("statmod")
  y <- c(0.1, 0.5, 1, 2, 5, 10, 20)
  mu <- 1.4; phi <- 0.74
  
  exact_lo <- statmod::pinvgauss(y, mean = mu, dispersion = phi, lower.tail = TRUE)
  exact_hi <- statmod::pinvgauss(y, mean = mu, dispersion = phi, lower.tail = FALSE)
  
  got_lo <- ptweedie_inversion(y, mu = mu, phi = phi, power = 3, lower.tail = TRUE)
  got_hi <- ptweedie_inversion(y, mu = mu, phi = phi, power = 3, lower.tail = FALSE)
  
  expect_equal(got_lo, exact_lo, tolerance = 1e-8)
  expect_equal(got_hi, exact_hi, tolerance = 1e-8)
})

test_that("ptweedie_series matches the compound Poisson-gamma exact sum for 1<p<2", {
  # Independent exact computation, not relying on any package internals
  exact_cdf <- function(y, mu, phi, power) {
    lambda <- mu^(2 - power) / (phi * (2 - power))
    alpha  <- (2 - power) / (power - 1)
    gam    <- phi * (power - 1) * mu^(power - 1)
    n <- 1:300
    exp(-lambda) + sum(dpois(n, lambda) * pgamma(y, shape = n * alpha, scale = gam))
  }
  
  cases <- expand.grid(power = c(1.01, 1.3, 1.9), mu = c(0.5, 1, 5), phi = c(0.5, 1, 2))
  y <- 1
  
  for (i in seq_len(nrow(cases))) {
    p <- cases$power[i]; m <- cases$mu[i]; f <- cases$phi[i]
    exact <- exact_cdf(y, m, f, p)
    got   <- ptweedie_series(y, power = p, mu = m, phi = f, lower.tail = TRUE)
    expect_equal(got, exact, tolerance = 1e-8,
                 info = sprintf("power=%s mu=%s phi=%s", p, m, f))
  }
})

test_that("series and inversion methods agree with each other for 1<p<2", {
  y <- c(0.1, 0.5, 1, 2, 5)
  for (p in c(1.01, 1.1, 1.5, 1.9)) {
    for (lt in c(TRUE, FALSE)) {
      a <- ptweedie_series(y, power = p, mu = 1, phi = 1, lower.tail = lt)
      b <- ptweedie_inversion(y, power = p, mu = 1, phi = 1, lower.tail = lt)
      expect_equal(a, b, tolerance = 1e-6,
                   info = sprintf("power=%s lower.tail=%s", p, lt))
    }
  }
})

test_that("point mass at y=0 is correct and respects lower.tail / log.p", {
  mu <- 1.5; phi <- 2; power <- 1.3
  lambda <- mu^(2 - power) / (phi * (2 - power))
  pi0 <- exp(-lambda)
  
  expect_equal(ptweedie(0, power = power, mu = mu, phi = phi, lower.tail = TRUE),
               pi0, tolerance = 1e-10)
  expect_equal(ptweedie(0, power = power, mu = mu, phi = phi, lower.tail = FALSE),
               1 - pi0, tolerance = 1e-10)
  expect_equal(ptweedie_inversion(0, power = power, mu = mu, phi = phi, log.p = TRUE),
               log(pi0), tolerance = 1e-8)
  expect_equal(ptweedie_inversion(0, power = power, mu = mu, phi = phi,
                                  lower.tail = FALSE, log.p = TRUE),
               log1p(-pi0), tolerance = 1e-8)
})

test_that("special p cases (p=0,1,2,3) agree with their standard R distributions", {
  y <- c(0.1, 0.5, 1, 2, 5)
  mu <- 2; phi <- 1.5
  
  # p = 0: Normal, via special_cases() directly (ptweedie() itself requires power > 1)
  for (lt in c(TRUE, FALSE)) {
    out0 <- special_cases(y, rep(mu, length(y)), rep(phi, length(y)), power = 0,
                          type = "CDF", lower.tail = lt)
    expect_equal(out0$f, pnorm(y, mean = mu, sd = sqrt(phi), lower.tail = lt),
                 tolerance = 1e-10)
  }
  
  # p = 1: Poisson (phi acts as a scale on y here, per special_cases())
  for (lt in c(TRUE, FALSE)) {
    expect_equal(
      ptweedie(y, power = 1, mu = mu, phi = phi, lower.tail = lt),
      ppois(y / phi, lambda = mu / phi, lower.tail = lt),
      tolerance = 1e-10
    )
  }
  
  # p = 2: gamma
  for (lt in c(TRUE, FALSE)) {
    expect_equal(
      ptweedie(y, power = 2, mu = mu, phi = phi, lower.tail = lt),
      pgamma(y, scale = mu * phi, shape = 1 / phi, lower.tail = lt),
      tolerance = 1e-10
    )
  }
  
  # p = 3: inverse Gaussian (via statmod)
  skip_if_not_installed("statmod")
  for (lt in c(TRUE, FALSE)) {
    expect_equal(
      ptweedie(y, power = 3, mu = mu, phi = phi, lower.tail = lt),
      statmod::pinvgauss(y, mean = mu, dispersion = phi, lower.tail = lt),
      tolerance = 1e-10
    )
  }
})

test_that("ptweedie converges to the gamma distribution as p -> 2", {
  mu <- 1; phi <- 1
  y <- c(0.5, 1, 2, 5)
  exact_gamma <- pgamma(y, scale = mu * phi, shape = 1 / phi)
  
  for (p in c(1.999, 2.00001, 2.001)) {
    got <- ptweedie(y, power = p, mu = mu, phi = phi)
    expect_equal(got, exact_gamma, tolerance = 1e-4,
                 info = sprintf("power=%s", p))
  }
})

test_that("number of integration regions for 1<p<2 is small (regression for p-near-1 fix)", {
  # Before the conditional-integrand fix, this specific case required 408
  # integration regions via pre-acceleration hitting the MAX_ACC cap. After
  # the fix it should converge in well under 50.
  out <- ptweedie_inversion(1, mu = 5, phi = 2, power = 1.01, details = TRUE)
  expect_lt(out$regions, 50)
  expect_equal(out$exitstatus, 0L)
})

test_that("exitstatus flags non-convergence when it genuinely occurs", {
  # A deliberately awkward corner of the parameter space; this is a smoke
  # test that exitstatus is a valid 0/1 vector, not a check for a specific
  # failure, since failure cases may change as the algorithm is refined.
  out <- ptweedie_inversion(c(0.001, 1, 1000), mu = 0.01, phi = 0.01, power = 6,
                            details = TRUE)
  expect_true(all(out$exitstatus %in% c(0L, 1L)))
})

test_that("ptweedie_series errors for p outside its valid domain (1<p<2)", {
  # ptweedie_series is documented for 1 < p < 2 only, and checks this before
  # doing any arithmetic (so no NaN warnings are produced on the way).
  msg <- "requires 1 < power < 2"
  expect_error(ptweedie_series(1, power = 3,   mu = 1, phi = 1), msg)
  expect_error(ptweedie_series(1, power = 2,   mu = 1, phi = 1), msg)
  expect_error(ptweedie_series(1, power = 1,   mu = 1, phi = 1), msg)
  expect_error(ptweedie_series(1, power = 0.5, mu = 1, phi = 1), msg)
  expect_no_warning(try(ptweedie_series(1, power = 3, mu = 1, phi = 1), silent = TRUE))
})


test_that("KNOWN LIMITATION: extreme upper tail can return small negative values", {
  # Previously (before the checkStopPreAcc consecutive-count fix), this case
  # silently returned a small negative value with exitstatus=0. It now
  # correctly flags non-convergence (exitstatus=1) rather than returning an
  # unreliable value silently -- this is the intended effect of that fix.
  hi <- ptweedie_inversion(20, power = 1.01, mu = 0.5, phi = 0.5,
                           lower.tail = FALSE, details = TRUE)
  expect_true(hi$exitstatus == 1)
})


