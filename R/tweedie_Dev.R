#' Tweedie deviance: Computes the unit deviance for Tweedie distributions
#'
#' Internal function to evaluate the Tweedie unit deviance.
#' \bold{Not intended for general users.}
#'
#' @param y vector of quantiles.
#' @param power The power parameter \eqn{p}{power}.
#' @param mu The mean parameter.
#' @return A numeric vector containing the unit deviance.
#' @keywords internal
#'
#' @noRd
tweedie_dev <- function(y, mu, power)
{
  # 
  # Peter K Dunn 
  # 29 Oct 2000 
  # 
  
  p <- power
  if(p == 1) {
    dev <- array( dim = length(y) )
    mu <- array( dim = length(y), mu )
    dev[y != 0] <- y[y != 0] * log(y[y != 0] / mu[y != 0]) - (y[y != 0] - mu[y != 0])
    dev[y == 0] <- mu[y == 0]
  } else{
    if(p == 2) {
      dev <- log(mu / y) + (y / mu) - 1
    } else{
      if (p == 0) {
        dev <- (y - mu) ^ 2
        dev <- dev / 2
      } else{
        dev <- (y ^ (2 - p)) / ((1 - p) * (2 - p)) - 
          (y * (mu ^ (1 - p)))/(1 - p) + 
          (mu ^ (2 - p)) / (2 - p)
      }
    }
  }
  dev * 2
}

