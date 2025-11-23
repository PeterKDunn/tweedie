
test_that("Deprecated functions work", {
  expect_no_error(
    ptweedie.inversion(0.001, mu = 0.01, phi = 0.01, power = 6, verbose=FALSE)  
  )
  expect_no_error(
    dtweedie.inversion(0.001, mu = 0.01, phi = 0.01, power = 6, verbose=FALSE)  
  )
})  
