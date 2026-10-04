# details = TRUE must work, and return full-length components, whether values
# come from special cases (exact formulas) or from the Fortran.

test_that("dtweedie_inversion details=TRUE works for special p (p = 3)", {
  out <- dtweedie_inversion(c(0.5, 1, 2), mu = 1, phi = 1, power = 3, details = TRUE)
  expect_type(out, "list")
  expect_length(out$exitstatus, 3)
  expect_equal(out$exitstatus, c(0L, 0L, 0L))
  expect_equal(out$density, statmod::dinvgauss(c(0.5, 1, 2), mean = 1, dispersion = 1))
})

test_that("dtweedie_inversion details=TRUE: exitstatus aligned when some y are special", {
  y <- c(0, 0.5, 1, 2)                       # y = 0 is a special case for 1 < p < 2
  out <- dtweedie_inversion(y, mu = 1, phi = 1, power = 1.5, details = TRUE)
  expect_length(out$exitstatus, length(y))
  expect_length(out$density, length(y))
  expect_equal(out$exitstatus[1], 0L)
  expect_equal(out$density, dtweedie_inversion(y, mu = 1, phi = 1, power = 1.5))
})

test_that("ptweedie_inversion details=TRUE works for special p (p = 3)", {
  out <- ptweedie_inversion(c(0.5, 1, 2), mu = 1, phi = 1, power = 3, details = TRUE)
  expect_length(out$exitstatus, 3)
  expect_equal(out$exitstatus, c(0L, 0L, 0L))
})

test_that("ptweedie_inversion details=TRUE: exitstatus aligned when some q are special", {
  q <- c(0, 0.5, 1, 2)
  out <- ptweedie_inversion(q, mu = 1, phi = 1, power = 1.5, details = TRUE)
  expect_length(out$exitstatus, length(q))
  expect_equal(out$exitstatus[1], 0L)
  expect_equal(out$cdf, ptweedie_inversion(q, mu = 1, phi = 1, power = 1.5))
})
