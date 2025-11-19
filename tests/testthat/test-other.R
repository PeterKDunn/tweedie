test_that("No common errors", {
  # These are situations that had casued errors in the past that 
  # I have (hopefully) now fixed
  expect_no_error(
    ptweedie.inversion(0.01, mu=5, phi=0.01, power=1.01, verbose=FALSE)  )
})

