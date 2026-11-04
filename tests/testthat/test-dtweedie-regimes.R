# dtweedie(): the regimes (since 3.2.0; no interpolation).
# Reference values: the series summed in Rmpfr (the Poisson-gamma series for
# 1 < p < 2, and the p > 2 series with enough bits to absorb its
# cancellation), or the inversion integral in 240-bit arithmetic; mu = 1.
# The cases: z = phi y^(p - 2) (the dispersion after rescaling y to 1) and
# y = exp(k sqrt(z)) for z < 1, exp(k) otherwise (k = -4, 4: both sides of
# the mean). The comment gives the method dtweedie() uses.

ref_cases <- read.table(header = FALSE, sep = ",", strip.white = TRUE,
                        col.names = c("p", "z", "k", "ref"), text = "
    1.01 , 1e-06, -4, 0.1329679670029987   # Corrected saddlepoint approximation
    1.01 , 1e-06,  4, 0.1346952652779497   # Corrected saddlepoint approximation
    1.01 , 0.01 , -4, 0.00062685630853999657   # Series
    1.01 , 0.01 ,  4, 0.0023165621284170564   # Series
    1.01 , 1    , -4, 1.0680004399812795e-19   # Series
    1.01 , 1    ,  4, 0.0012095343102612316   # Series
    1.01 , 100  , -4, 1.3176141328965586e-155   # Series
    1.01 , 100  ,  4, 6.9297316227483992e-159   # Series
    1.5  , 1e-06, -4, 0.13436624783131737   # Corrected saddlepoint approximation
    1.5  , 1e-06,  4, 0.1332956061255457   # Corrected saddlepoint approximation
    1.5  , 0.01 , -4, 0.0019420826276924685   # Corrected saddlepoint approximation
    1.5  , 0.01 ,  4, 0.00087263397533110705   # Corrected saddlepoint approximation
    1.5  , 1    , -4, 0.00031047785341921532   # Series
    1.5  , 1    ,  4, 1.0415371661354074e-07   # Series
    1.5  , 100  , -4, 0.018791815584854925   # Series
    1.5  , 100  ,  4, 6.303951839154811e-06   # Series
    1.99 , 1e-06, -4, 0.13577714916571051   # Corrected saddlepoint approximation
    1.99 , 1e-06,  4, 0.131908466917949   # Corrected saddlepoint approximation
    1.99 , 0.01 , -4, 0.0051556038086133932   # Corrected saddlepoint approximation
    1.99 , 0.01 ,  4, 0.00028166469458795296   # Corrected saddlepoint approximation
    1.99 , 1    , -4, 0.91300205691668979   # Series
    1.99 , 1    ,  4, 9.0722139749698276e-24   # Series
    1.99 , 100  , -4, 0.48827443622553968   # Series
    1.99 , 100  ,  4, 0.00010447782088521249   # Series
    2.001, 1e-06, -4, 0.13580896801295506   # Corrected saddlepoint approximation
    2.001, 1e-06,  4, 0.13187747003734346   # Corrected saddlepoint approximation
    2.001, 0.01 , -4, 0.0052614615972436402   # Corrected saddlepoint approximation
    2.001, 0.01 ,  4, 0.0002740645187974158   # Corrected saddlepoint approximation
    2.001, 1    , -4, 0.98890671170832556   # Inversion
    2.001, 1    ,  4, 1.6601613596138186e-24   # Inversion
    2.001, 100  , -4, 0.50528322740740739   # Inversion
    2.001, 100  ,  4, 0.00010621174287158586   # Inversion
    2.05 , 1e-06, -4, 0.13595078437489844   # Corrected saddlepoint approximation
    2.05 , 1e-06,  4, 0.13173946909170689   # Corrected saddlepoint approximation
    2.05 , 0.01 , -4, 0.0057554362703466349   # Corrected saddlepoint approximation
    2.05 , 0.01 ,  4, 0.00024235451305556319   # Corrected saddlepoint approximation
    2.05 , 1    , -4, 1.3725668892336513   # Inversion
    2.05 , 1    ,  4, 3.8420088514243133e-28   # Inversion
    2.05 , 100  , -4, 0.58293554452706975   # Series
    2.05 , 100  ,  4, 0.00011232553324736351   # Series
    2.5  , 1e-06, -4, 0.13725914137086068   # Corrected saddlepoint approximation
    2.5  , 1e-06,  4, 0.13047790863812864   # Corrected saddlepoint approximation
    2.5  , 0.01 , -4, 0.012382173530526426   # Corrected saddlepoint approximation
    2.5  , 0.01 ,  4, 7.1612513534260382e-05   # Corrected saddlepoint approximation
    2.5  , 1    , -4, 7.2182876085070298   # Series
    2.5  , 1    ,  4, 7.5979881142752412e-114   # Series
    2.5  , 100  , -4, 1.3552345024489267   # Series
    2.5  , 100  ,  4, 3.5696431051314019e-05   # Series
    4    , 1e-06, -4, 0.14169876223667979   # Corrected saddlepoint approximation
    4    , 1e-06,  4, 0.12634744797337605   # Corrected saddlepoint approximation
    4    , 0.01 , -4, 0.085883908156664787   # Inversion
    4    , 0.01 ,  4, 2.8382226886549565e-07   # Inversion
    4    , 1    , -4, 19.827331994772798   # Series
    4    , 100  , -4, 3.48731806962659   # Series
    4    , 100  ,  4, 8.5017882921583924e-233   # Series
    8    , 1e-06, -4, 0.15414520835965076   # Corrected saddlepoint approximation
    8    , 1e-06,  4, 0.11587857517322507   # Corrected saddlepoint approximation
    8    , 0.01 , -4, 1.0620576947207279   # Inversion
    8    , 0.01 ,  4, 1.4513522695624838e-23   # Inversion
    8    , 1    , -4, 28.303845460930585   # Series
    8    , 100  , -4, 6.8906381514942714   # Series
    12   , 1e-06, -4, 0.16751704649830215   # Corrected saddlepoint approximation
    12   , 1e-06,  4, 0.10616621211135716   # Corrected saddlepoint approximation
    12   , 0.01 , -4, 2.6617588903990437   # Series
    12   , 0.01 ,  4, 3.8429248657301771e-85   # Series
    12   , 1    , -4, 33.308540296100126   # Series
    12   , 100  , -4, 9.3190400077408313  # Series
")

test_that("dtweedie() matches high-precision reference values in every regime", {
  y   <- exp(ref_cases$k * pmin(sqrt(ref_cases$z), 1))
  phi <- ref_cases$z / y^(ref_cases$p - 2)
  got <- mapply(function(y, phi, p) dtweedie(y, mu = 1, phi = phi, power = p), 
                y, phi, ref_cases$p)
  relerr <- abs(got / ref_cases$ref - 1)
  expect_lt(max(relerr), 1e-12)
})

test_that("dtweedie() gives the same values for a vector as one at a time", {
  set.seed(1)
  for (p in c(1.05, 1.5, 2.5, 8)) {
    y   <- c(0, exp(rnorm(30, sd = 3)))
    mu  <- exp(rnorm(31))
    phi <- exp(rnorm(31, sd = 2))
    f   <- dtweedie(y, mu = mu, phi = phi, power = p)
    f1  <- mapply(function(y, mu, phi) dtweedie(y, mu = mu, phi = phi, power = p), y, mu, phi)
    expect_equal(f, f1, tolerance = 1e-14)
  }
})

test_that("dtweedie() is accurate at small dispersion (where it once interpolated)", {
  # Before 3.2.0 the interpolation had relative errors up to about 1e-5 here
  # (p = 1.01, z = 1e-8) and about 2.5e-7 for 7 <= p <= 10. 
  # Reference: the inversion integral in 240-bit arithmetic.
  phi <- 1e-8
  expect_lt(abs(dtweedie(1, mu = 1, phi = phi, power = 1.01) / 3989.4228006733515031 - 1), 1e-13)
})

test_that("dtweedie() handles y = 0, y = Inf, NA, and names", {
  expect_equal(dtweedie(0, mu = 2, phi = 1, power = 1.5), exp(-2^0.5 / 0.5))
  expect_equal(dtweedie(0, mu = 2, phi = 1, power = 2.5), 0)
  expect_equal(dtweedie(c(Inf, NA), mu = 1, phi = 1, power = 1.5), c(0, NA))
  expect_equal(dtweedie(c(Inf, NA), mu = 1, phi = 1, power = 2.5), c(0, NA))
  expect_equal(names(dtweedie(c(a = 1, b = 2), mu = 1, phi = 1, power = 1.5)), c("a", "b"))
  expect_equal(dtweedie(numeric(0), mu = 1, phi = 1, power = 1.5), numeric(0))
})

test_that("the corrected saddlepoint approximation reduces to Stirling's series at p = 2, and is exact at p = 3", {
  # p = 2 (gamma, shape a = 1/z): f / f_sp = sqrt(2 pi) a^(a - 1/2) exp(-a) / Gamma(a)
  for (z in c(0.01, 0.05)) {
    a     <- 1 / z
    exact <- exp(0.5 * log(2 * pi) + (a - 0.5) * log(a) - a - lgamma(a))
    sp    <- tweedie:::dtweedie_saddle_corrected(1, 1, z, 2 + 1e-12)
    fsp   <- (2 * pi * z)^(-1/2)
    expect_lt(abs(sp$density / fsp / exact - 1), 1e-12)
  }
  # p = 3: the correction is zero (the saddlepoint approximation is exact)
  sp <- tweedie:::dtweedie_saddle_corrected(c(0.5, 1, 2), 1, 0.01, 3)
  expect_equal(sp$density, statmod::dinvgauss(c(0.5, 1, 2), mean = 1, dispersion = 0.01), tolerance = 1e-14)
})

test_that("tweedie_dev() has full relative accuracy near y = mu", {
  # For p = 2 the deviance is 2 {(y/mu - 1) - log(y/mu)}, computed here from
  # its series (with y/mu - 1 = e exactly: e is a power of 2)
  for (e in 2^c(-10, -20, -30)) {
    y   <- 1 + e
    ref <- 2 * sum( (-1)^(2:30) * e^(2:30) / (2:30) )       # = 2 {e - log1p(e)}
    expect_lt(abs(tweedie_dev(y, mu = 1, power = 2) / ref - 1), 1e-14)
  }
  # For p = 3/2 the deviance is 4 (sqrt(y) - sqrt(mu))^2, written here without
  # cancellation as 4 {(y - mu)/(sqrt(y) + sqrt(mu))}^2
  for (y in c(1 + 2^-27, 1 - 2^-20, 2, 1e-3)) {
    ref <- 4 * ( (y - 1) / (sqrt(y) + 1) )^2
    expect_lt(abs(tweedie_dev(y, mu = 1, power = 1.5) / ref - 1), 1e-14)
  }
})

test_that("dtweedie_inversion() uses the method it is given", {
  y  <- c(0.5, 2)
  d1 <- dtweedie_inversion(y, mu = 1, phi = 1, power = 2.5, method = 1, details = TRUE)
  d3 <- dtweedie_inversion(y, mu = 1, phi = 1, power = 2.5, method = 3, details = TRUE)
  expect_equal(d1$method, c(1, 1))
  expect_equal(d3$method, c(3, 3))
  expect_equal(d1$density, d3$density, tolerance = 1e-10)
  expect_error(dtweedie_inversion(1, mu = 1, phi = 1, power = 2.5, method = 4))
})

test_that("dtweedie_series() does not need a huge matrix for widely spread y (1 < p < 2)", {
  # Before 3.2.0 one range of N served all values of y, and this needed about 15 Gb
  y <- c(1e-10, 1e6)
  f <- dtweedie_series(y, power = 1.5, mu = 1, phi = 0.01)
  expect_true(all(is.finite(f)))
  expect_equal(as.vector(f), c(dtweedie_series(1e-10, power = 1.5, mu = 1, phi = 0.01),
                               dtweedie_series(1e6,   power = 1.5, mu = 1, phi = 0.01)))
})
