test_that("No common errors", {
  # These are situations that had caused errors in the past that 
  # I have (hopefully) now fixed
  expect_no_error(
    dtweedie_inversion(power=1.5, mu=1, phi=0.5513086, y=0.5557088)
  )
})

