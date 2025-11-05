test_that("Negative y-values return 0 for density", {
  fy <- dtweedie(-1, mu = 1, phi = 1, power = 3.5)
  expect_equal(fy, 0)
})


test_that("Negative y-values return 0 for DF", {
  fy <- ptweedie(-1, mu = 1, phi = 1, power = 3.5)
  expect_equal(fy, 0)

  fy <- ptweedie(-1, mu = 1, phi = 1, power = 1.5)
  expect_equal(fy, 0)
})




test_that("Invalid parameters cause error in dtweedie", {
  # 1. Test for bad mu value (mu <= 0)
  # The error message is the string "mu must be positive."
  expect_error(
    dtweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  # Assuming the error message is "phi must be positive."
  expect_error(
    dtweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
})




test_that("Invalid parameters cause error in ptweedie", {
  
  # 1. Test for bad mu value (mu <= 0)
  # The error message is the string "mu must be positive."
  expect_error(
    ptweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  # Assuming the error message is "phi must be positive."
  expect_error(
    ptweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
})


test_that("Invalid parameters cause error in rtweedie", {
  
  # 1. Test for bad mu value (mu <= 0)
  # The error message is the string "mu must be positive."
  expect_error(
    rtweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  # Assuming the error message is "phi must be positive."
  expect_error(
    rtweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
})


test_that("Invalid parameters cause error in qtweedie", {
  
  # 1. Test for bad mu value (mu <= 0)
  # The error message is the string "mu must be positive."
  expect_error(
    qtweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  # Assuming the error message is "phi must be positive."
  expect_error(
    qtweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
})

