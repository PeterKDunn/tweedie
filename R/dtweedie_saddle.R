#' Tweedie densities evaluation using the saddlepoint approximation
#' 
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
#' @importFrom base pi
#'
#' @aliases dtweedie.saddle
#' 
#' @references
#' 	Nelder, J. A. and Pregibon, D. (1987).
#' 	An extended quasi-likelihood function
#' 	\emph{Biometrika},
#' 	\bold{74}(2), 221--232.
#' 	\doi{10.1093/biomet/74.2.221}

#' 
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


#' @export
dtweedie.saddle <- function(y, xi = NULL, mu, phi, eps = 1/6, power = NULL){ 
  .Deprecated("dtweedie_saddle", package = "tweedie")
  dtweedie_saddle(y, xi = NULL, mu, phi, eps = 1/6, power = NULL)
}


