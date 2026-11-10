#' @title Tweedie Distribution: Convert Between Parameter Formats
#' @name tweedie_convert
#' @description
#' For \eqn{1 < p < 2}, a Tweedie random variable is a Poisson sum of gamma
#' random variables: \eqn{Y = \sum_{i=1}^N Z_i}{Y = Z_1 + ... + Z_N}, where
#' \eqn{N \sim \mbox{Pois}(\lambda)}{N ~ Pois(lambda)} and the \eqn{Z_i} are
#' independent gamma variables with shape \eqn{\alpha}{alpha} and scale
#' \eqn{\gamma}{gamma}.
#'
#' \code{tweedie_convert()} converts from the Tweedie (\acronym{glm})
#' parameters \eqn{p}, \eqn{\mu}{mu} and \eqn{\phi}{phi} to the underlying
#' Poisson parameter \eqn{\lambda}{lambda} and gamma parameters
#' \eqn{\alpha}{alpha} and \eqn{\gamma}{gamma}.
#'
#' \code{tweedie_convert_back()} converts the other way: from
#' \eqn{\lambda}{lambda}, \eqn{\alpha}{alpha} and \eqn{\gamma}{gamma}
#' to \eqn{p}, \eqn{\mu}{mu} and \eqn{\phi}{phi}.
#'
#' @details
#' The relationships are
#' \deqn{\lambda = \frac{\mu^{2 - p}}{\phi (2 - p)}, \quad
#'       \alpha  = \frac{2 - p}{p - 1}, \quad
#'       \gamma  = \phi (p - 1) \mu^{p - 1}}{lambda = mu^(2 - p) / (phi * (2 - p)),
#'       alpha = (2 - p) / (p - 1),
#'       gamma = phi * (p - 1) * mu^(p - 1)}
#' and, in reverse,
#' \deqn{p = \frac{\alpha + 2}{\alpha + 1}, \quad
#'       \mu = \lambda \alpha \gamma, \quad
#'       \phi = \frac{\lambda \alpha (\alpha + 1) \gamma^2}{\mu^p}.}{p = (alpha + 2) / (alpha + 1),
#'       mu = lambda * alpha * gamma,
#'       phi = lambda * alpha * (alpha + 1) * gamma^2 / mu^p.}
#'
#' Here \eqn{\gamma}{gamma} is the gamma \emph{scale} parameter; the
#' corresponding \emph{rate} parameter (sometimes written \eqn{\beta}{beta})
#' is \eqn{\beta = 1/\gamma}{beta = 1/gamma}.
#' In \R's notation, the gamma variables are
#' \code{rgamma(n, shape = alpha, scale = gamma)}.
#'
#' @param xi a synonym for \code{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param power the power parameter \eqn{p}{power}; a synonym for \eqn{\xi}{xi}.
#' @param lambda the mean \eqn{\lambda}{lambda} of the underlying Poisson distribution.
#' @param alpha the shape parameter \eqn{\alpha}{alpha} of the underlying gamma distribution.
#' @param gamma the scale parameter \eqn{\gamma}{gamma} of the underlying gamma distribution
#'   (the reciprocal of the rate parameter).
#'
#' @return
#' \code{tweedie_convert()} returns a list of the parameters of the
#' corresponding underlying Poisson and gamma distributions:
#' \code{poisson.lambda} (\eqn{\lambda}{lambda} from the underlying Poisson distribution),
#' \code{gamma.shape}, \code{gamma.scale} (the shape \eqn{\alpha}{alpha} and
#'  scale \eqn{\gamma}{gamma} parameters from the underlying gamma distribution),
#'  \code{p0} (the probability that \eqn{Y = 0}),
#'  \code{gamma.mean} and \code{gamma.phi} (the gamma mean \eqn{\alpha\gamma}{alpha * gamma}
#'  and dispersion parameter \eqn{1/\alpha}{1/alpha} values).
#'
#' \code{tweedie_convert_back()} returns a list containing
#' \code{mu} (the Tweedie mean \eqn{\mu}{mu}),
#' \code{phi} (the Tweedie dispersion parameter \eqn{\phi}{phi}),
#' \code{power} (the Tweedie index parameter \eqn{p}),
#' and \code{p0} (the probability that \eqn{Y = 0}).
#'
#' @importFrom stats glm fitted
#' @importFrom statmod tweedie
#'
#' @examples
#' ### Fit a Tweedie density
#' pretend <- data.frame( y = rgamma(20, shape = 1, rate = 1) )
#' fit <- glm(y ~ 1, data = pretend,
#'            family = statmod::tweedie(link.power = 0, var.power = 1.4))
#'
#' # Convert parameters
#' tweedie_convert(mu = fitted(fit, type="response"), phi = 1, power = 1.4)
#'
#' ### Convert, then convert back
#' pg <- tweedie_convert(mu = 2, phi = 0.5, power = 1.4)
#' tweedie_convert_back(lambda = pg$poisson.lambda,
#'                      alpha  = pg$gamma.shape,
#'                      gamma  = pg$gamma.scale)
#'
#' @export
tweedie_convert <- function(xi = NULL, mu, phi, power = NULL){
  ### ADDED 14 July 2017

  # SORT OUT THE NOTATION (i.e., xi VS power)
  out            <- sort_notation(xi = xi, power = power)
  xi             <- out$xi
  power          <- out$power
  index.par.long <- out$index.par.long

  # Error checks
  if ( !is.numeric(power) ) stop( paste(index.par.long, "must be numeric.\n") )
  if ( power <= 1)     stop( paste(index.par.long, "must be greater than 1.\n") )
  if ( power >= 2)     stop( paste(index.par.long, "must be less than 2.\n") )
  if ( any(phi <= 0) ) stop("phi must be positive.\n")
  if ( any(mu <= 0) )  stop("mu must be positive.\n")

  if ( length(mu) != length(phi) ){
    if ( length(mu)  == 1 ) mu  <- rep_len(mu,  length(phi) )
    if ( length(phi) == 1 ) phi <- rep_len(phi, length(mu) )
  }
  # Now  mu  and  phi  will be the same length if one of them happened to be a scalar, so this works:
  if ( length(mu) != length(phi) ) stop("phi and mu must be scalars, or the same length.\n")

  lambda <- ( mu ^ (2 - power) ) / ( phi * (2 - power) )  # Poisson distribution mean
  alpha  <- (2 - power)  / (power - 1)                    # gamma distribution shape
  gam    <- phi * (power - 1) * mu ^ (power - 1)          # gamma distribution scale
  p0     <- exp( -lambda )
  mu.g   <- alpha * gam                                   # gamma distribution mean
  phi.g  <- 1 / alpha                                     # gamma distribution dispersion

  list( poisson.lambda = lambda,
        gamma.shape    = alpha,
        gamma.scale    = gam,
        p0             = p0,
        gamma.mean     = mu.g,
        gamma.phi      = phi.g)
}



#' @rdname tweedie_convert
#' @export
tweedie_convert_back <- function(lambda, alpha, gamma){

  # Error checks
  if ( !is.numeric(lambda) | !is.numeric(alpha) | !is.numeric(gamma) ) {
    stop("lambda, alpha and gamma must be numeric.\n")
  }
  if ( any(lambda <= 0) ) stop("lambda must be positive.\n")
  if ( any(alpha  <= 0) ) stop("alpha must be positive.\n")
  if ( any(gamma  <= 0) ) stop("gamma must be positive.\n")

  # Inputs must be scalars, or all the same length
  len <- max( length(lambda), length(alpha), length(gamma) )
  if ( !all( c(length(lambda), length(alpha), length(gamma)) %in% c(1, len) ) ) {
    stop("lambda, alpha and gamma must be scalars, or the same length.\n")
  }
  lambda <- rep_len(lambda, len)
  alpha  <- rep_len(alpha,  len)
  gamma  <- rep_len(gamma,  len)

  power <- (alpha + 2) / (alpha + 1)                           # Tweedie index parameter
  mu    <- lambda * alpha * gamma                              # Tweedie mean
  phi   <- lambda * alpha * (alpha + 1) * gamma ^ 2 / mu ^ power   # Tweedie dispersion
  p0    <- exp( -lambda )

  # If alpha is a scalar, return power as a scalar too
  if ( length(unique(power)) == 1 ) power <- power[1]

  list( mu    = mu,
        phi   = phi,
        power = power,
        p0    = p0)
}



#' @rdname tweedie_convert
#' @export
tweedie.convert <- function(xi = NULL, mu, phi, power = NULL){
  lifecycle::deprecate_warn(when = "3.0.5",
                            what = "tweedie.convert()",
                            with = "tweedie_convert()")
  if (is.null(power)) power <- xi
  tweedie_convert(xi    = NULL,
                  mu    = mu,
                  phi   = phi,
                  power = power)
}
