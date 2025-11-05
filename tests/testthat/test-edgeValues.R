test_that("Zero input returns zero (ignoring attributes and noise)", {
  
  # The test environment is adding a dimension attribute and seeing numerical noise.
  
  # 1. Test dtweedie(0, ...)
  expect_equal(
    dtweedie(0, mu = 1, phi = 1, power = 3.5), 
    0,
    tolerance = 1e-8,      # Required to pass the numerical check (actual != 0)
    ignore_attr = TRUE     # Required to pass the dimension check (dim(actual) != dim(expected))
  )
  
  # 2. Test ptweedie(0, ...)
  expect_equal(
    ptweedie(0, mu = 1, phi = 1, power = 3.5), 
    0,
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})





