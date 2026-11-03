# When the pre-acceleration zone reaches its cap of 400 regions, acceleration
# is started before tmax, where it is not valid, and can appear to converge
# to a wrong value. Since 3.2.0 this is flagged (exitstatus = 1).

test_that("ptweedie_inversion() flags a result when pre-acceleration hits its cap", {
  # F underflows to 0 here (log F is about -5e12), but the inversion returns 1.2e-8
  expect_warning(o <- ptweedie_inversion(0.001, mu = 0.01, phi = 0.01, power = 6,
                                         details = TRUE),
                 "did not reach the target accuracy")
  expect_equal(o$exitstatus, 1L)
  expect_gt(o$regions, 400)
})

test_that("ptweedie() is unaffected: it gets these cases right without a warning", {
  # p > 2, xi = 2.5e9: the saddlepoint approximation is used
  expect_silent(F <- ptweedie(0.001, mu = 0.01, phi = 0.01, power = 6))
  expect_equal(F, 0)
  # p > 2, xi <= 1e5 (p = 4, mu = 0.1, phi = 0.01: xi = 2500): the flagged
  # inversion value is tiny, so it is recomputed by integrating the density
  expect_silent(F <- ptweedie(0.001, mu = 0.1, phi = 0.01, power = 4))
  expect_equal(F, 0)
  # 1 < p < 2: the series is used
  expect_silent(F <- ptweedie(1, mu = 15, phi = 0.01, power = 1.01))
  expect_equal(F, 0)
})

test_that("ptweedie(): a flagged inversion value near one is replaced by one minus the other tail", {
  # p = 2.5, mu = 0.01, phi = 0.001 (xi = 2e4), y = mu/2: the inversion hits
  # the pre-acceleration cap; S(y) is 1 to double precision. This used to
  # warn after the cap was flagged; it should not.
  expect_silent(S <- ptweedie(5e-3, mu = 0.01, phi = 0.001, power = 2.5, lower.tail = FALSE))
  expect_equal(S, 1)
  expect_silent(lS <- ptweedie(5e-3, mu = 0.01, phi = 0.001, power = 2.5, lower.tail = FALSE, log.p = TRUE))
  expect_lte(lS, 0)
  expect_gt(lS, -1e-300)
})
