# The regimes (R/tweedie_regimes.R): one place for the thresholds, and
# functions that say which method each value goes to.

test_that("tweedie_regime_settings() holds the thresholds", {
  S <- tweedie:::tweedie_regime_settings()
  expect_equal(S$d_tol_saddle, 1e-13)
  expect_equal(S$d_tol_series, 1e-12)
  expect_equal(c(S$p_xi_trust, S$p_xi_body_int, S$p_xi_body_sp, S$p_xi_saddle), 
               c(1e4, 1e5, 1e6, 1e7))
})

test_that("dtweedie_regime() gives the method dtweedie() uses", {
  # p = 1.5: z = phi y^(p - 2); the corrected saddlepoint for z below about 0.07
  r <- tweedie:::dtweedie_regime(c(0, 1, 1, NA, Inf), mu = 1, phi = c(1, 1e-3, 1, 1, 1), power = 1.5)
  expect_equal(r, c("exact", "saddlepoint", "series", NA, "exact"))
  r <- tweedie:::dtweedie_regime(c(1, 1), mu = 1, phi = c(1e-6, 1), power = 4)
  expect_equal(r, c("saddlepoint", "series/inversion"))
  expect_equal(tweedie:::dtweedie_regime(1:3, mu = 1, phi = 1, power = 3), rep("exact", 3))
  # and dtweedie() uses the corrected saddlepoint exactly there
  y   <- exp(seq(-3, 3, length = 41))
  r   <- tweedie:::dtweedie_regime(y, mu = 1, phi = 0.05, power = 2.5)
  sp  <- tweedie:::dtweedie_saddle_corrected(y, mu = 1, phi = 0.05, power = 2.5)$density
  f   <- dtweedie(y, mu = 1, phi = 0.05, power = 2.5)
  expect_equal(f[r == "saddlepoint"], sp[r == "saddlepoint"])
  expect_true(any(r == "saddlepoint") && any(r != "saddlepoint"))
})

test_that("ptweedie_regime() follows xi = mu^(2-p)/{phi (p-2)}", {
  p   <- 2.5
  xi  <- c(1, 1e4, 1e4 * 1.01, 1e5 * 1.01, 1e7 * 1.01)
  phi <- 1 / (xi * (p - 2))
  r   <- tweedie:::ptweedie_regime(rep(1, 5), mu = 1, phi = phi, power = p)
  expect_equal(r, c("inversion", "inversion", "integration in tails", "integration", "saddlepoint"))
  # at xi = 1e7: the saddlepoint when phi mu^(p - 2) <= 1e-5 (p = 2.5), not near p = 2
  expect_equal(tweedie:::ptweedie_regime(1, mu = 1, phi = 1 / (1e7 * 0.5), power = 2.5), "saddlepoint")
  expect_equal(tweedie:::ptweedie_regime(1, mu = 1, phi = 1 / (1e7 * 0.001), power = 2.001), "integration")
  expect_equal(tweedie:::ptweedie_regime(c(0, 1), mu = 1, phi = 1, power = 1.5), c("exact", "series"))
})
