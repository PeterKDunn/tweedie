#' Density for the Tweedie Family of Distributions using Saddlepoint Approximation
#'
#' Evaluates the probability density function (PDF) for Tweedie distributions
#' for given values of the dependent variable \code{y}, the mean \code{mu},
#' dispersion \code{phi}, and power parameter \code{power}.
#'
#' @description The Tweedie family of distributions belong to the class of
#'   exponential dispersion models (EDMs), famous for their role in generalized
#'   linear models. The Tweedie distributions are the EDMs with a variance
#'   of the form \eqn{\mbox{var}[Y] = \phi\mu^p}{var[Y] = phi*mu^p} where \eqn{p \ge 1}{p >= 1}.
#'   \bold{This function only evaluates for \eqn{p \ge 1}{p >= 1}.}
#'
#' @details
#' Evaluation is difficult for \eqn{p}{p} outside of \eqn{p=0, 1, 2, 3}{power = 0, 1, 2, 3}. This function
#' uses one of two primary methods, depending on the combination of parameters:
#' \enumerate{
#'   \item Evaluation of an infinite series (\code{dtweedie_series}).
#'   \item Interpolation from stored values computed via a Fourier inversion technique (\code{dtweedie_inversion}).
#' }
#' This function employs a two-dimensional interpolation procedure to compute
#' the density for some parts of the parameter space from previously computed
#' values (interpolation) and uses the series solution for others.
#'
#' Special cases include the Poisson (\eqn{p = 1} with \eqn{\phi = 1}{phi = 1}), gamma (\eqn{p = 2}), and inverse Gaussian (\eqn{p = 3}) distributions.
#'
#' @section Note:
#' \code{dtweedie} and \code{ptweedie} are the only functions meant for users. Consequently, all checks on the function inputs are performed here.
#'
#' @param y Vector of quantiles.
#' @param power A synonym for \eqn{\xi}{xi}, the Tweedie index parameter.
#' @param mu The mean parameter.
#' @param phi The dispersion parameter.
#'
#' @return The density, approximated using the saddlepoint approximation.
#'
#' @seealso \code{\link{ptweedie}}, \code{\link{rtweedie}}, \code{\link{dtweedie_series}}, \code{\link{dtweedie_inversion}}, \code{\link{dtweedie_saddle}}
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
#' Jorgensen, B. (1997).
#' \emph{Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' @aliases dtweedie_stable
#' @aliases dtweedie.stable
#'
#' @export
dtweedie_stable <- function(y, power, mu, phi)
{
  # Error checks
  if ( power < 1) stop("power must be greater than 2.\n")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y < 0) ) stop("y must be a non-negative vector.\n")
  if ( any(mu <= 0) ) stop("mu must be positive.\n")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.\n")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.\n")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  density <- y
  alpha <- (2 - power) / (1 - power)
  beta <- 1
  k <- 1  # The parameterization used
  delta <- 0
  gamma <- phi * (power - 1) * 
    ( 1 / (phi * (power - 2)) * cos( alpha * pi / 2 ) ) ^ (1 / alpha)
  
  #        require(stabledist) # Needed for  dstable
  ds <- stabledist::dstable(y,
                            alpha = alpha, 
                            beta = beta, 
                            gamma = gamma, 
                            delta = delta, 
                            pm = k)
  density <- exp((y * mu ^ (1 - power) / (1 - power) - mu ^ (2 - power) / (2 - power)) / phi) * ds
  
  
  density
}



#############################################################################
