test_that("Negative y-values return 0 for density", {
  expect_equal(
      dtweedie(-1, mu = 1, phi = 1, power = 3.5),
      0
  )
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

  # 3. Test for bad n value (n > 0)
  expect_error(
    rtweedie(0, mu = 1, phi = -1, power = 1.5),
    regexp = "n must be a positive integer."
  )
})


test_that("Invalid parameters cause error in qtweedie", {
  
  # 1. Test for bad mu value (mu <= 0)
  expect_error(
    qtweedie(1, mu = -1, phi = 1, power = 1.5), 
    regexp = "mu must be positive." 
  )
  expect_error(
    qtweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  expect_error(
    qtweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
  expect_error(
    qtweedie(1, mu = 1, phi = -1, power = 3.5),
    regexp = "phi must be positive."
  )
  
  # 3. Test for q less than zero
  expect_error(
    qtweedie(-1, mu = 1, phi = -1, power = 1.5),
    regexp = "quantiles must be between 0 and 1."
  )
  expect_error(
    qtweedie(-1, mu = 1, phi = -1, power = 3.5),
    regexp = "quantiles must be between 0 and 1."
  )
  
  
  # 4. Test for q greater than zero
  expect_error(
    qtweedie(2, mu = 1, phi = -1, power = 1.5),
    regexp = "quantiles must be between 0 and 1."
  )
  expect_error(
    qtweedie(2, mu = 1, phi = -1, power = 3.5),
    regexp = "quantiles must be between 0 and 1."
  )
  
  
  # 5. Test for q greater than zero, less than 0  
  expect_error(
    qtweedie(c(-1, 2), mu = 1, phi = -1, power = 1.5),
    regexp = "quantiles must be between 0 and 1."
  )
  expect_error(
    qtweedie(c(-1, 2), mu = 1, phi = -1, power = 3.5),
    regexp = "quantiles must be between 0 and 1."
  )
})


