
test_that("Deprecated functions work", {
  expect_warning(
    ptweedie.inversion(1, mu = 1, phi = 1, power = 6, verbose=FALSE)  
  )
  expect_warning(
    dtweedie.inversion(1, mu = 1, phi = 1, power = 6, verbose=FALSE)  
  )
  expect_warning(
    dtweedie.series(1, mu = 1, phi = 1, power = 6)  
  )
  expect_warning(
    ptweedie.series(1, mu = 1, phi = 1, power = 1.5)  
  )
})  
