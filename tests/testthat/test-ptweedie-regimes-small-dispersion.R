# ptweedie() for p > 2: the methods used at small dispersion (since 3.2.0).
# Reference values: the inversion integral evaluated in 240-bit arithmetic
# (as in the paper), with mu = 1 and xi = mu^(2-p) / {phi (p-2)}.

test_that("ptweedie() is accurate in the body for 1e5 < xi < 1e7 near p = 2 (density integrated)", {
  # The inversion's relative error in the body reaches about 1e-8 here; the
  # integration's is about 1e-10
  p <- 2.001; phi <- 1 / (1e6 * (p - 2))
  expect_equal(ptweedie(1, mu = 1, phi = phi, power = p),
               0.5042073462989560273, tolerance = 1e-10)
  expect_equal(ptweedie(0.936754446796629, mu = 1, phi = phi, power = p),
               0.021015626545158896249, tolerance = 1e-10)
  expect_equal(ptweedie(1.06324555320337, mu = 1, phi = phi, power = p, lower.tail = FALSE),
               0.024430716672937429629, tolerance = 1e-10)
  phi <- 1 / (3e6 * (p - 2))
  expect_equal(ptweedie(0.95801795735350004, mu = 1, phi = phi, power = p),
               0.010000000264980043324, tolerance = 1e-10)
})

test_that("ptweedie() uses the saddlepoint approximation at xi = 1e7 when phi mu^(p-2) <= 1e-5", {
  # p = 2.01: phi mu^(p-2) = 1e-5. The inversion's error here was about 1e-8
  p <- 2.01; phi <- 1 / (1e7 * (p - 2))
  y <- c(0.99367544467966296, 1, 1.0063245553203399)
  expect_equal(ptweedie(y, mu = 1, phi = phi, power = p),
               tweedie:::ptweedie_saddle(y, mu = 1, phi = phi, power = p))
  expect_equal(ptweedie(1, mu = 1, phi = phi, power = p),
               0.50042262471581144823, tolerance = 1e-9)
  expect_equal(ptweedie(1.0063245553203399, mu = 1, phi = phi, power = p, lower.tail = FALSE),
               0.022921447176722724304, tolerance = 1e-9)
  # p = 2.001: phi mu^(p-2) = 1e-4, so the saddlepoint (error about 2e-8
  # here) is not used at xi = 1e7; the tail is integrated
  p <- 2.001; phi <- 1 / (1e7 * (p - 2))
  S  <- ptweedie(1.08, mu = 1, phi = phi, power = p, lower.tail = FALSE)
  expect_equal(S, 3.1122261997262424e-15, tolerance = 1e-10)
  expect_gt(abs(S / tweedie:::ptweedie_saddle(1.08, mu = 1, phi = phi, power = p, lower.tail = FALSE) - 1), 1e-9)
})

test_that("ptweedie() replaces a failed inversion near p = 2 by a rough value, with a warning", {
  # p = 2.05, xi = 1 (phi = 20), very small y: the inversion fails (it used to
  # return 0.5, or 0, with a warning); the saddlepoint approximation is rough
  # but of the right size. The value is still flagged.
  p <- 2.05; phi <- 1 / (1 * (p - 2))
  expect_warning(F <- ptweedie(1.35924654511923e-36, mu = 1, phi = phi, power = p),
                 "did not reach the target accuracy")
  expect_lt(abs(F / 1.2467929812465447816e-26 - 1), 0.1)
  expect_warning(F <- ptweedie(1e-10, mu = 1, phi = phi, power = p),
                 "did not reach the target accuracy")
  expect_true(F > 0.05 && F < 0.2)                    # the true value is 0.0924
})
