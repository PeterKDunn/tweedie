test_that("No common errors", {
  # These are situations that had caused errors in the past that 
  # I have (hopefully) now fixed
  expect_no_error(
    ptweedie_inversion(0.01, mu=5, phi=0.01, power=1.01, verbose=FALSE)  )
  expect_no_error(
    ptweedie_inversion(0.005005, mu=1, phi=1, power=3, verbose=FALSE)  )
  expect_no_error(
    ptweedie_inversion(0.095095, mu=1, phi=1, power=4.5, verbose=FALSE)  )
  expect_no_error(
    ptweedie_inversion(0.001, mu=.01, phi=0.01, power=6, verbose=FALSE)  )
})

