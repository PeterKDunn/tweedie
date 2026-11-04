#' @title Tweedie densities evaluation using the saddlepoint approximation
#' @name dtweedie_saddle
#' @description
#' Density function for the Tweedie EMDs using a saddlepoint approximation.
#' 
#' @param y vector of quantiles.
#' @param xi scalar; the value of \eqn{\xi}{xi} such that the variance is \eqn{\mbox{var}[Y]=\phi\mu^{\xi}}{var[Y] = phi * mu^xi}. A synonym for \code{power}.
#' @param mu vector of mean \eqn{\mu}{mu}.
#' @param phi vector of dispersion parameters \eqn{\phi}{phi}.
#' @param power scalar; a synonym for \eqn{\xi}{xi}, the Tweedie index parameter.
#' @param eps the offset in computing the variance function; the default is \code{eps=1/6} (as suggested by Nelder and Pregibon, 1987).
#' 
#' @return A numeric vector of densities.
#' 
#' @examples
#' # Plot a Tweedie density
#' y <- seq(0.01, 4, length = 50)
#' fy <- dtweedie_saddle(y, power = 1.1, mu = 1, phi = 1)
#' plot(y, fy, type = "l", lwd = 2, ylab = "Density")
#'
#'
#' @keywords distribution
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#' 
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' 	Nelder, J. A. and Pregibon, D. (1987).
#' 	An extended quasi-likelihood function
#' 	\emph{Biometrika},
#' 	\bold{74}(2), 221--232.
#' 	\doi{10.1093/biomet/74.2.221}
#'
#' @keywords distribution
#' @export
dtweedie_saddle <- function(y, xi = NULL, mu, phi, eps = 1/6, power = NULL) {

  # Sort out the xi/power notation
  if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
  xi.notation <- TRUE
  if ( is.null(power) ) {
    power <- xi
  } else {
    xi.notation <- FALSE
  }
  if ( is.null(xi) ) {
    xi.notation <- FALSE
    xi <- power
  }
  if ( xi != power ) {
    cat("Different values for xi and power given; the value of xi used.\n")
    power <- xi
  }
  index.par       <- ifelse( xi.notation, "xi", "p")
  index.par.long  <- ifelse( xi.notation, "xi", "power")
  
  
  # Error traps
  #
  if( any(phi <= 0) )
    stop("phi must be positive.")
  if( (power >= 1) & any(y < 0))
    stop("y must be a non-negative vector.")
  if( any(mu <= 0) )
    stop("mu must be positive.")
  if ( length(mu) == 1 )
    mu <- array( dim=length(y), mu )
  if ( length(phi) == 1 )
    phi <- array( dim=length(y), phi )
  #
  # END error checks
  #
  #   if ( any(y==0) ) y <- y + (1/6)
  y.eps <- y 
  
  if (power<2) y.eps <- y + eps
  # Nelder and Pregibon's idea for problems at y=0
  
  y0 <- (y == 0)
  density <- y
  #    if(any(y == 0)) {
  #        if(power >= 2) {
  #         # This is a problem!
  #            density[y0] <- 0 * y[y0]
  #        }
  #      else{
  
  dev <- tweedie_dev(y = y, 
                     mu = mu, 
                     power = power)
  density <- (2 * pi * phi * y.eps ^ power) ^ (-1/2) * exp( -dev / (2 * phi) )
  
  #      }
  #   }
  if ( any(y == 0) ){
    if((power >= 1) && (power < 2)) {
      lambda <- mu[y0] ^ (2 - power) / (phi[y0] * (2 - power))
      density[y0] <- exp( -lambda)
    } else {
      density[y0] <- 0
    }
  }
  
  #    if(any(y != 0)) {
  #     dev <- tweedie_dev(y=y[yp], mu=mu[yp], power=power)
  #      density[yp] <- (2*pi*phi[yp]*y[yp]^power)^(-1/2) * exp( -dev/(2*phi[yp]) )
  #    }
  density
}


#' @rdname dtweedie_saddle
#' @export
dtweedie.saddle <- function(y, xi = NULL, mu, phi, eps = 1/6, power = NULL){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "dtweedie.saddle()", 
                            with = "dtweedie_saddle()")
  dtweedie_saddle(y, xi = NULL, mu, phi, eps = 1/6, power = NULL)
}



#' @noRd
dtweedie_saddle_corrected <- function(y, mu, phi, power) {
  # The saddlepoint approximation to the density, with its asymptotic
  # correction (used by dtweedie() where the dispersion is small; see there).
  # For y > 0 (and p not 0, 1, 2 or 3, where dtweedie() uses closed forms).
  #
  # The saddlepoint approximation is
  #   f_sp(y) = (2 pi phi y^p)^(-1/2) exp{-d(y, mu)/(2 phi)},
  # with d the unit deviance. Its relative error depends only on p and
  #   z = phi y^(p - 2),
  # (the dispersion after rescaling y to 1), and has the asymptotic expansion
  #   f(y) / f_sp(y) ~ 1 + a_1(p) z + a_2(p) z^2 + ...
  # (the Edgeworth expansion of the exponentially tilted density at its mean;
  # e.g. Barndorff-Nielsen and Cox, 1989). Here a_1(p) = p (p - 3)/24, and each
  # a_m(p) is a polynomial in p of degree 2m. They were derived in exact
  # rational arithmetic from the cumulants of the tilted distribution, whose
  # standardised cumulants are rho_r = z^(r/2 - 1) prod_{j=1}^{r-2} {j (p - 1) + 1};
  # as checks, a_m(3) = 0 for every m (the saddlepoint approximation is exact
  # for the inverse Gaussian), and a_m(2) are the coefficients of the series
  # for Gamma(a) / {sqrt(2 pi) a^(a - 1/2) exp(-a)} (Stirling's series; gamma
  # distribution, z = 1/a): -1/12, 1/288, 139/51840, -571/2488320, ...
  # Each a_m(p) is symmetric about p = 3/2 (a_m(p) = a_m(3 - p)), and is stored
  # below as a polynomial in s = (p - 3/2)^2 (which also evaluates it without
  # cancellation; in powers of p it would not be).
  #
  # The series diverges, so it is summed only while its terms decrease (at most
  # eight terms); the size of the last term summed is returned as an estimate
  # of the relative error ('relerr'), which is reliable: against high-precision
  # reference values, the actual error is below the estimate (or below about
  # 3e-14, whichever is larger).
  #
  # For 1 < p < 2 the density also has a term that the expansion cannot
  # capture: Y is a Poisson sum of gamma variables, and when the number of
  # gamma terms near the peak, N = 1/{z (2 - p)}, is small compared with their
  # shape, alpha = (2 - p)/(p - 1), the density has ripples (one for each N).
  # Their relative size is about exp(-2 pi^2 N / alpha), which is added to the
  # error estimate (it matters only for p near 1).
  #
  # The unit deviance is computed without cancellation (tweedie_dev()), and
  # the density on the log scale, so it neither cancels nor underflows early.
  
  n   <- max(length(y), length(mu), length(phi))
  y   <- rep_len(y,   n)
  mu  <- rep_len(mu,  n)
  phi <- rep_len(phi, n)
  
  sc   <- saddle_correction(z = phi * y ^ (power - 2), 
                            power = power)
  corr   <- sc$corr
  relerr <- sc$relerr
  logf <- -0.5 * log(2 * pi * phi) - (power / 2) * log(y) -
            tweedie_dev(y = y, mu = mu, power = power) / (2 * phi) + 
            log(pmax(corr, 0))
  
  list(density = exp(logf),
       logf    = logf,
       relerr  = relerr)
}


#' @noRd
saddle_correction <- function(z, power) {
  # The asymptotic correction to the saddlepoint approximation to the density,
  # f / f_sp ~ 1 + a_1(p) z + a_2(p) z^2 + ..., with z = phi y^(p - 2), and an
  # estimate of the relative error of the corrected approximation (see
  # dtweedie_saddle_corrected()). It depends only on z and p, so it also
  # decides, before any density is computed, where dtweedie() uses the
  # corrected saddlepoint approximation (see dtweedie_regime()).
  
  coefs <- list(                    # a_m(p): coefficients of s^0, s^1, ..., s^m
    c(-0.09375, 0.041666666666666664),
    c(-0.00732421875, 0.048177083333333329, -0.019965277777777776),
    c(-0.0016021728515625, 0.034427897135416663, -0.075949435763888878, 0.027095389660493824),
    c(-0.00056326389312744141, 0.02581825256347656, -0.14382620917426214, 0.19916336624710645, -0.062351581589184669),
    c(-0.00027107074856758118, 0.021643268636294771, -0.22847773044828384, 0.70417957793678132, -0.72854336273803344, 0.20391846525419138),
    c(-0.00016518373740836978, 0.020416948771370309, -0.35073252484143252, 1.9103887863270608, -4.0759231478858871, 3.414915257654259, -0.86700313188385614),
    c(-0.00012167552085884381, 0.021555557473870324, -0.54841883637046152, 4.6489037289541795, -16.734359918660804, 27.459283327147066, -19.506648265746609, 4.542651684388753),
    c(-0.00010504019574142376, 0.025265580179458472, -0.8931267961325523, 10.869733400438863, -58.99721957176417, 157.9847748331855, -211.72895167670441, 131.35384807074067, -28.305925780171091)
  )
  
  n  <- length(z)
  s  <- (power - 3 / 2) ^ 2
  am <- sapply(coefs, function(cf) {
    v <- 0
    for (k in length(cf):1) v <- v * s + cf[k]
    v
  })
  
  # The correction: sum the terms while they decrease in size
  corr   <- rep(1, n)
  last   <- rep(Inf, n)
  going  <- rep(TRUE, n)
  zm     <- rep(1, n)
  for (m in seq_along(am)) {
    zm     <- zm * z
    term   <- am[m] * zm
    going  <- going & (abs(term) < last)
    corr[going] <- corr[going] + term[going]
    last[going] <- abs(term[going])
  }
  relerr <- last
  relerr[last == Inf] <- 0          # (z = 0 cannot occur for y > 0; for safety)
  
  if ( (power > 1) && (power < 2) ) {
    N_peak <- 1 / ( z * (2 - power) )
    alpha  <- (2 - power) / (power - 1)
    relerr <- relerr + exp(-2 * pi^2 * N_peak / alpha)
  }
  
  # (A correction factor that is not positive means z is far too large for
  # the expansion: no value)
  bad <- !is.na(corr) & (corr <= 0)
  relerr[bad] <- Inf
  list(corr   = corr,
       relerr = relerr)
}
