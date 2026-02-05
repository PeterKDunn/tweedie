test_that("Works OK", {
  # For upstream dependencies
  expect_no_error({
    y <- rgamma(5, 1, 1)
    m1 <- glm(y ~ 1, 
              family = statmod::tweedie(var.power = 2.5, 
                                        link.power = 0))
    logLiktweedie(m1)}
  )}
)

