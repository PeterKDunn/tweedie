# For 1 < p < 2, ptweedie() must give small tail probabilities to full relative
# accuracy, by using the series rather than the inversion (whose error is
# absolute, about 1e-15). References: 50-60 digit sums of the series (mpmath).

test_that("ptweedie: small upper tails for 1 < p < 2 are relatively accurate", {
  y   <- c(50, 100, 200, 300)
  ref <- c(2.234650455582427e-10, 3.4136489462303752e-20,
           2.4362963524509583e-40, 8.2746323426280853e-61)
  v <- ptweedie(y, mu = 1, phi = 4, power = 1.5, lower.tail = FALSE)
  expect_lt(max(abs(v - ref) / ref), 1e-12)
  lv <- ptweedie(y, mu = 1, phi = 4, power = 1.5, lower.tail = FALSE, log.p = TRUE)
  expect_lt(max(abs(lv - log(ref))), 1e-11)
})

test_that("ptweedie: small lower tails for 1 < p < 2 are relatively accurate", {
  ref <- 1.0093541776661044e-19
  v <- ptweedie(0.3, mu = 1, phi = 0.01, power = 1.5)
  expect_lt(abs(v - ref) / ref, 1e-12)
})

test_that("ptweedie: moderate probabilities for 1 < p < 2 are unchanged", {
  y <- seq(0.05, 10, length = 50)
  expect_lt(max(abs(ptweedie(y, mu = 1, phi = 4, power = 1.5) - pchisq(y, df = 0, ncp = 1))), 1e-14)
})
