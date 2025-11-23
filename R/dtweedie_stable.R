#' Tweedie Distribution: Stable Evaluation for the Density
#'
#' Internal function to evaluate the Tweedie density using stable distributions.
#' \bold{Not intended for general users.}
#'
#' @param y vector of quantiles.
#' @param power The power parameter \eqn{p}{power}.
#' @param mu The mean parameter.
#' @param phi The dispersion parameter.
#' @return A numeric vector of densities.
#' @keywords internal

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
