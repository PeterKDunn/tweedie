# For p < 4/3, Im k(t) oscillates about zero before the -t y term takes over,
# so a region between the zeros found can hold many oscillations of the
# integrand. A single Gauss rule over such a region was wrong, without a
# warning (errors up to 2e-3 before 3.2.1); the regions are now integrated
# adaptively. Also, when the tail regions are only rounding noise, the
# acceleration now stops (it used to run to 400 regions and warn, on some
# platforms only).

test_that("inversion, p near 1, small y: agrees with the series", {
  cases <- data.frame(
    y   = c(1e-4, 1e-4, 1e-4, 1e-3, 1e-4, 1e-3, 1e-4, 1e-3, 1e-4, 1e-4, 1e-4),
    mu  = c(1,    1,    5,    5,    0.01, 1,    5,    15,   15,   5,    15),
    phi = c(0.1,  1,    1,    10,   10,   10,   10,   10,   10,   10,   10),
    p   = c(1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.05, 1.2,  1.2))
  for (i in seq_len(nrow(cases))) with(cases[i, ], {
    ref <- ptweedie_series(y, mu = mu, phi = phi, power = p)
    out <- ptweedie_inversion(y, mu = mu, phi = phi, power = p, details = TRUE)
    expect_equal(out$exitstatus, 0L)
    # (the inversion promises absolute accuracy)
    expect_lt(abs(out$cdf - ref), 1e-13)
    S <- ptweedie_inversion(y, mu = mu, phi = phi, power = p, lower.tail = FALSE)
    expect_lt(abs(S - (1 - ref)), 1e-13)
  })
})

test_that("inversion converges, without a warning, when the tail is negligible", {
  # F(0.001) here is the point mass at zero, exp(-lambda)
  expect_no_warning(
    F <- ptweedie_inversion(0.001, mu = 5, phi = 10, power = 1.01)
  )
  expect_equal(F, exp(-5^0.99 / (10 * 0.99)), tolerance = 1e-14)
})
