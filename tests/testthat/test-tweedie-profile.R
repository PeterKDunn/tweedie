# Tests for tweedie_profile()
#
# The data sets are small, and only a few values of p are used, so these run
# quickly. tweedie_profile() prints a progress bar, so calls are wrapped in
# quiet_profile() (or capture_output()) to keep the test output clean.

quiet_profile <- function(...) {
  out <- NULL
  utils::capture.output( out <- suppressWarnings( tweedie_profile(...) ) )
  out
}

# Data with exact zeros (so 1 < p < 2) and without (so p > 1 is allowed)
set.seed(2026)
n.sim  <- 150
x.sim  <- runif(n.sim)
mu.sim <- exp(1 + x.sim)
dat.zeros <- data.frame(x = x.sim,
                        y = rtweedie(n.sim, mu = mu.sim, phi = 2, power = 1.6))
dat.pos   <- data.frame(x = x.sim,
                        y = rtweedie(n.sim, mu = mu.sim, phi = 0.3, power = 2.6))

data(rock)



test_that("The output has the expected components, in either notation", {
  out.xi <- quiet_profile(perm ~ 1, data = rock,
                          xi.vec = seq(1.5, 2.5, length = 5))
  expect_true( all(c("x", "y", "ht", "L", "xi", "xi.max", "L.max",
                     "phi", "phi.max", "ci", "method", "phi.method") %in% names(out.xi)) )
  expect_length(out.xi$L, 5)
  expect_length(out.xi$phi, 5)
  expect_length(out.xi$ci, 2)

  # Supplying  p.vec  (not  xi.vec) returns  p  and  p.max  instead
  out.p <- quiet_profile(perm ~ 1, data = rock,
                         p.vec = seq(1.5, 2.5, length = 5))
  expect_true( all(c("p", "p.max") %in% names(out.p)) )
  expect_false( any(c("xi", "xi.max") %in% names(out.p)) )
  expect_equal(out.p$p.max, out.xi$xi.max)
})



test_that("Without smoothing, the MLEs are the best of the values computed", {
  xi.vec <- seq(1.3, 1.9, by = 0.1)
  out <- quiet_profile(y ~ x, data = dat.zeros, xi.vec = xi.vec,
                       do.smooth = FALSE)
  best <- which.max(out$L)
  expect_equal(out$xi.max,  xi.vec[best])
  expect_equal(out$L.max,   max(out$L))
  expect_equal(out$phi.max, out$phi[best])
  expect_true( all(is.finite(out$L)) )
  expect_true( all(out$phi > 0) )
})



test_that("Each L is the log-likelihood at that p and phi", {
  # Refit each GLM from scratch (no warm start), and recompute the
  # log-likelihood directly with dtweedie(). This checks both that L
  # (reused from the search for phi) is the log-likelihood, and that starting
  # each glm.fit() from the previous fit gives the same fitted values.
  # (p = 2 checks the gamma special case.)
  xi.vec <- c(1.5, 1.75, 2, 2.25, 2.5)
  out <- quiet_profile(perm ~ 1, data = rock, xi.vec = xi.vec,
                       do.smooth = FALSE)
  for (i in seq_along(xi.vec)) {
    fit <- glm(perm ~ 1, data = rock,
               family = statmod::tweedie(var.power = xi.vec[i], link.power = 0))
    logL <- sum( dtweedie(rock$perm, mu = fitted(fit), phi = out$phi[i],
                          power = xi.vec[i], log = TRUE) )
    expect_equal(out$L[i], logL, tolerance = 1e-6)
  }
})



test_that("Each phi maximises the log-likelihood at its p", {
  xi.vec <- c(1.4, 1.6, 1.8)
  out <- quiet_profile(y ~ x, data = dat.zeros, xi.vec = xi.vec,
                       do.smooth = FALSE)
  for (i in seq_along(xi.vec)) {
    mu <- fitted( glm(y ~ x, data = dat.zeros,
                      family = statmod::tweedie(var.power = xi.vec[i],
                                                link.power = 0)) )
    logL <- function(phi) sum( dtweedie(dat.zeros$y, mu = mu, phi = phi,
                                        power = xi.vec[i], log = TRUE) )
    # Moving phi away from the estimate (either way) lowers the log-likelihood
    expect_gt(out$L[i], logL(out$phi[i] * 0.98))
    expect_gt(out$L[i], logL(out$phi[i] * 1.02))
  }
})



test_that("The profile finds p close to the true value, inside the CI", {
  # Data simulated with p = 1.6 (dat.zeros) and p = 2.6 (dat.pos)
  out1 <- quiet_profile(y ~ x, data = dat.zeros,
                        xi.vec = seq(1.2, 1.9, by = 0.1))
  expect_true( abs(out1$xi.max - 1.6) < 0.2 )
  expect_true( out1$ci[1] < out1$xi.max & out1$xi.max < out1$ci[2] )

  out2 <- quiet_profile(y ~ x, data = dat.pos,
                        xi.vec = seq(2.2, 3.0, by = 0.2))
  expect_true( abs(out2$xi.max - 2.6) < 0.3 )
  expect_true( out2$phi.max > 0 )
})



test_that("method only matters for the saddlepoint approximation", {
  # With phi.method = "mle", L comes from dtweedie() while finding phi,
  # so the  series  and  inversion  methods give the same profile
  xi.vec <- seq(1.3, 1.9, by = 0.2)
  out.inv <- quiet_profile(y ~ x, data = dat.zeros, xi.vec = xi.vec,
                           do.smooth = FALSE, method = "inversion")
  out.ser <- quiet_profile(y ~ x, data = dat.zeros, xi.vec = xi.vec,
                           do.smooth = FALSE, method = "series")
  expect_equal(out.ser$L,   out.inv$L)
  expect_equal(out.ser$phi, out.inv$phi)

  # The saddlepoint approximation gives a (slightly) different profile
  out.sad <- quiet_profile(y ~ x, data = dat.pos, xi.vec = c(2.2, 2.6, 3.0),
                           do.smooth = FALSE, method = "saddlepoint")
  expect_equal(out.sad$phi.method, "saddlepoint")
  expect_true( all(is.finite(out.sad$L)) )
})



test_that("phi.method = 'saddlepoint' gives the mean deviance", {
  xi.vec <- c(1.5, 2, 2.5)
  out <- quiet_profile(perm ~ 1, data = rock, xi.vec = xi.vec,
                       do.smooth = FALSE, phi.method = "saddlepoint")
  for (i in seq_along(xi.vec)) {
    mu <- fitted( glm(perm ~ 1, data = rock,
                      family = statmod::tweedie(var.power = xi.vec[i],
                                                link.power = 0)) )
    phi.dev <- mean( tweedie_dev(y = rock$perm, mu = mu, power = xi.vec[i]) )
    expect_equal(out$phi[i], phi.dev, tolerance = 1e-6)
  }
})



test_that("Impossible values of p are removed, with a message", {
  # p <= 1 (except p = 0) is never possible
  expect_output(
    out <- suppressWarnings(
      tweedie_profile(perm ~ 1, data = rock, do.smooth = FALSE,
                      xi.vec = c(0.5, 1.5, 2, 2.5))),
    regexp = "have been removed"
  )
  expect_equal(out$xi, c(1.5, 2, 2.5))

  # With exact zeros in the response, p must be between 1 and 2
  expect_output(
    out <- suppressWarnings(
      tweedie_profile(y ~ x, data = dat.zeros, do.smooth = FALSE,
                      xi.vec = c(1.4, 1.6, 2.2, 2.6))),
    regexp = "must be between 1 and 2"
  )
  expect_equal(out$xi, c(1.4, 1.6))
})



test_that("Default values of p depend on whether the response has zeros", {
  out.zeros <- quiet_profile(y ~ 1, data = dat.zeros, do.smooth = FALSE)
  expect_equal(out.zeros$xi, seq(1.2, 1.8, by = 0.1))

  out.pos <- quiet_profile(perm ~ 1, data = rock, do.smooth = FALSE)
  expect_equal(out.pos$xi, seq(1.5, 5, by = 0.5))
})



test_that("Bad inputs give errors or warnings", {
  expect_error(
    capture_output(
      tweedie_profile(perm ~ 1, data = rock, xi.vec = c(1.5, 2),
                      do.smooth = FALSE, conf.level = 1.5)),
    regexp = "Confidence level"
  )

  # Smoothing needs at least five values of p
  expect_warning(
    capture_output(
      tweedie_profile(perm ~ 1, data = rock, xi.vec = c(1.5, 2, 2.5))),
    regexp = "at least five"
  )
})



test_that("fit.glm = TRUE returns the GLM fitted at the MLE of p", {
  out <- quiet_profile(perm ~ area, data = rock,
                       xi.vec = c(1.5, 1.75, 2, 2.25, 2.5),
                       do.smooth = FALSE, fit.glm = TRUE)
  expect_true( "glm.obj" %in% names(out) )
  # (Use the same convergence criterion as tweedie_profile(): IRLS converges
  # slowly here, so the default  epsilon  gives slightly different estimates)
  fit <- glm(perm ~ area, data = rock,
             control = glm.control(epsilon = 1e-09),
             family = statmod::tweedie(var.power = out$xi.max, link.power = 0))
  expect_equal(unname(out$glm.obj$coefficients), unname(coef(fit)),
               tolerance = 1e-6)
})



test_that("Offsets and add0 are accepted", {
  dat <- dat.pos
  dat$expo <- runif(nrow(dat), 1, 2)
  out <- quiet_profile(y ~ x + offset(log(expo)), data = dat,
                       xi.vec = c(2.2, 2.6, 3.0), do.smooth = FALSE)
  expect_true( all(is.finite(out$L)) )

  # add0 = TRUE also profiles p = 0 (the normal distribution)
  out0 <- quiet_profile(perm ~ 1, data = rock,
                        xi.vec = c(1.5, 2, 2.5), do.smooth = FALSE,
                        add0 = TRUE)
  expect_equal(out0$xi[1], 0)
  expect_length(out0$L, 4)
})



test_that("Prior weights give each observation dispersion phi/weight", {
  # With every weight equal to 2, the MLE of phi doubles (phi/2 is unchanged),
  # and the profile log-likelihood is unchanged
  xi.vec <- seq(1.5, 2.5, by = 0.25)
  rock.w <- rock
  rock.w$w2 <- 2
  out1 <- quiet_profile(perm ~ 1, data = rock.w, xi.vec = xi.vec,
                        do.smooth = FALSE)
  capture_output(
    out2 <- suppressWarnings(
      tweedie_profile(perm ~ 1, data = rock.w, weights = w2, xi.vec = xi.vec,
                      do.smooth = FALSE)) )
  expect_equal(out2$L, out1$L, tolerance = 1e-6)
  expect_equal(out2$phi / 2, out1$phi, tolerance = 1e-4)

  # An observation with zero weight is the same as leaving it out
  rock.w$w0 <- c(0, rep(1, nrow(rock) - 1))
  capture_output(
    out3 <- suppressWarnings(
      tweedie_profile(perm ~ 1, data = rock.w, weights = w0, xi.vec = xi.vec,
                      do.smooth = FALSE)) )
  out4 <- quiet_profile(perm ~ 1, data = rock[-1, ], xi.vec = xi.vec,
                        do.smooth = FALSE)
  expect_equal(out3$L, out4$L)
})



test_that("With smoothing, the MLE is the maximum of the spline", {
  out <- quiet_profile(perm ~ 1, data = rock,
                       xi.vec = seq(1.5, 2.75, length = 11))
  # The MLE need not be one of the values of p computed (or of the 50 grid
  # points), and is at least as high as all of them
  expect_true( out$L.max >= max(out$y) )
  expect_true( out$L.max >= max(out$L) )
  expect_true( out$ci[1] < out$xi.max & out$xi.max < out$ci[2] )
  # The CI limits are where the profile is  ht
  ss <- stats::splinefun(out$xi, out$L)
  expect_equal(ss(out$ci), rep(out$ht, 2), tolerance = 1e-4)
})



test_that("Inputs are checked", {
  # A misspelt method is an error (it used to silently use "series")
  expect_error(
    capture_output(
      tweedie_profile(perm ~ 1, data = rock, xi.vec = c(1.5, 2),
                      do.smooth = FALSE, method = "inverson")) )

  # Supplying both  p.vec  and  xi.vec  gives a warning
  expect_warning(
    capture_output(
      tweedie_profile(perm ~ 1, data = rock, xi.vec = c(1.5, 2),
                      p.vec = c(1.6, 2.1), do.smooth = FALSE)),
    regexp = "only  xi.vec  is used"
  )

  # No possible values of p is an error
  expect_error(
    capture_output(
      tweedie_profile(perm ~ 1, data = rock, xi.vec = c(0.5, 1),
                      do.smooth = FALSE)),
    regexp = "No values"
  )
})



test_that("The output is consistent with the values of p actually used", {
  # Impossible values are removed, so  p  must match  L  and  phi
  capture_output(
    out <- suppressWarnings(
      tweedie_profile(perm ~ 1, data = rock, p.vec = c(0.5, 1, 1.5, 2, 2.5),
                      do.smooth = FALSE)) )
  expect_equal(out$p, c(1.5, 2, 2.5))
  expect_length(out$L, 3)
  expect_length(out$phi, 3)

  # A single value of p works
  out1 <- quiet_profile(perm ~ 1, data = rock, xi.vec = 2, do.smooth = FALSE)
  expect_equal(out1$xi.max, 2)

  # Plotting does not change the results
  grDevices::pdf(NULL)
  out.plot <- quiet_profile(perm ~ 1, data = rock, xi.vec = c(1.5, 2, 2.5),
                            do.smooth = FALSE, do.plot = TRUE)
  grDevices::dev.off()
  out.noplot <- quiet_profile(perm ~ 1, data = rock, xi.vec = c(1.5, 2, 2.5),
                              do.smooth = FALSE)
  expect_equal(out.plot, out.noplot)
})
