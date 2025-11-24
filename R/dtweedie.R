#' Density for the Tweedie Family of Distributions
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
#' @param xi The value of \eqn{\xi}{xi} such that the variance is \eqn{\mbox{var}[Y]=\phi\mu^{\xi}}{var[Y] = phi * mu^xi}. A synonym for \code{power}.
#' @param mu The mean parameter.
#' @param phi The dispersion parameter.
#' @param power A synonym for \eqn{\xi}{xi}, the Tweedie index parameter.
#' @param verbose Logical. Whether to display details of the internal process. Defaults to \code{FALSE}.
#' @param details Logical. Whether to return computational reports (regions needed, etc.). Defaults to \code{FALSE}.
#'
#' @return
#' If \code{details = FALSE} (default), a numeric vector of densities.
#' If \code{details = TRUE}, a list containing:
#' \describe{
#'   \item{\code{density}}{The computed density vector.}
#'   \item{\code{regions}}{The number of integration regions used.}
#' }
#'
#' @importFrom stats dgamma dnorm dpois
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
#' @examples
#' # The examples from your .Rd file go here.
#' ### Plot a Tweedie density
#' power <- 2.5
#' mu <- 1
#' phi <- 1
#' y <- seq(0, 6, length = 500)
#' fy <- dtweedie(y = y, power = power, mu = mu, phi = phi)
#' plot(y, fy, type = "l", lwd = 2, ylab = "Density")
#' # Compare to the saddlepoint density
#' f.saddle <- dtweedie_saddle( y = y, power = power, mu = mu, phi = phi)
#' lines( y, f.saddle, col = 2 )
#' legend("topright", col = c(1,2), lwd = c(2,1),
#'   legend = c("Actual", "Saddlepoint") )
#'
#' @export
dtweedie <- function(y, xi = NULL, mu, phi, power = NULL, verbose = FALSE, details = FALSE){

  # Evaluates the density for Tweedie distributions, for given values of:
  #   y (possibly a vector)
  #   mu, the mean 
  #   phi the dispersion parameter
  #   power,  the Tweedie index parameter
  #   verbose: whether to display what is happening
  #   details:  whether to returns reports of relerr, regions needed, etc.

  # Two methods are employed:  cgf inversion (type=1)
  # and series evaluation (type=2 if 1<p<2; type=3 if p>2)).
  #
  # This function uses bivariate interpolation to accurately
  # approximate the inversion in conjunction with the series.
  #
  #
  # FIRST, establish whether the series or the cgf inversion is the better method.
  
  #
  # Here is a summary of what happens:
  #
  #   p                          |  What happpens
  # -----------------------------+-------------------
  #  1 -> p.lowest (=1.3)        |  Use series A
  #  p.lowest (=1.3) -> p=2      |  Use series A if
  #                              |        xix > xix.smallp (=0.8)
  #                              |      inversion with rho if
  #                              |        xix < xix.smallp (=0.8)
  #  p=2 -> p=3                  |  Use series B if
  #                              |        xix > xix.smallp (=0.8)
  #                              |      inversion with rho if
  #                              |        xix < xix.smallp (=0.8)
  #  p=3 -> p.highest (=4)       |  Use series B if xix > xix.bigp (=0.9)
  #                              |      inversion with 1/rho if
  #                              |        xix < xix.bigp (=0.9)
  #  p > 4                       |  Use series B if xix > 0.03
  #                              | Saddlepoint approximation if xix<= 0.03 (NEEDS FIXING)
  #
  
  ### BEGIN preliminary work
  
  # SORT OUT THE NOTATION (i.e., xi VS power)
  if (verbose) cat("- Checking notation\n")
  out <- sort_notation(xi = xi, power = power)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long ### MAY NOT BE NEEDED!!!
  
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(y, mu, phi, power)
  mu <- out$mu
  phi <- out$phi

  # density2 is the whole vector.
  # density is just the part where !special_y_cases.
  # All is resolved in the end
  density2 <- numeric(length = length(q) )
  if (details) regions <- integer(length = length(q))
  
  
  # IDENTIFY SPECIAL CASES
  special_y_cases <- rep(FALSE, length(y))
  if (verbose) cat("- Checking for special cases\n")
  out <- special_cases(y, mu, phi, power,
                       type = "PDF")

  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases

  if (verbose & special_p_cases) cat("  - Special case for p used\n")
  if ( any(special_y_cases) ) {
    special_y_cases <- out$special_y_cases  
    if (verbose) cat("  - Special cases for first input found\n")
#    density[special_y_cases] <- out$f[special_y_cases]
    density2 <- out$f
  }
  
  # Set things up for the interpolation/series; i.e., not special_y_cases
  density <- density2[ !special_y_cases ]
  mu <- mu[ !special_y_cases ]
  phi <- phi[ !special_y_cases ]
  y <- y[ !special_y_cases ]
  ### END preliminary work
  

  regions = NA  # 'regions'  may not be relevant if interpolation is used 
  if (special_p_cases){
    density = out$f
  } else {
    # Set up
    id.type0 <- array( FALSE, dim = length(y) )
    id.series <- id.type0
    id.interp <- id.type0
    density <- density2[ !special_y_cases ]
    
    ###   Now consider the cases 1 < p < 2   ###
    id.type0 <- (y == 0)
    if (any(id.type0) ) {
      if (power > 2) {
        density[id.type0] <- 0
      } else {
        lambda <- mu[id.type0] ^ (2 - power) / (phi[id.type0] * (2 - power))
        density[id.type0] <- exp( -lambda )
      }
    }
    
    
    xi <- array( dim = length(y) )
    xi[id.type0] <- 0
    xi[!id.type0] <- phi[!id.type0] * y[!id.type0] ^ (power - 2)
    xix <- xi / ( 1 + xi )

    if ( (power > 1) && (power <= 1.1) ) {
      id.series <- (!id.type0)
      if (any(id.series)){
        density[id.series] <- dtweedie_series(y = y[id.series],
                                              mu = mu[id.series], 
                                              phi = phi[id.series],
                                              power = power)
      }
      return(density = density)
    }
    
    if ( power == 1 ) { # AND phi not equal to one here
      id.series <- rep(TRUE, length(id.series))
      id.interp <- rep(FALSE, length(id.series))
    }
    if ( (power > 1.1) && (power <= 1.2) ) {
      id.interp <- ( (xix > 0) & (xix < 0.1) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 1.1
        p.hi <- 1.2
        xix.lo <- 0
        xix.hi <- 0.1
        np <- 15
        nx <- 25
      }
    }
    
    if ( (power > 1.2) && (power <= 1.3) ) {
      id.interp <- ( (xix > 0) & (xix < 0.3) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 1.2
        p.hi <- 1.3
        xix.lo <- 0
        xix.hi <- 0.3
        np <- 15
        nx <- 25
      }
    }
    if ( (power > 1.3) && (power <= 1.4) ) {
      id.interp <- ( (xix > 0) & (xix < 0.5) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 1.3
        p.hi <- 1.4
        xix.lo <- 0
        xix.hi <- 0.5
        np <- 15
        nx <- 25
      }
    }
    if ( (power > 1.4) && (power <= 1.5) ) {
      id.interp <- ( (xix > 0) & (xix < 0.8) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 1.4
        p.hi <- 1.5
        xix.lo <- 0
        xix.hi <- 0.8
        np <- 15
        nx <- 25
      }
    }
    if ( (power > 1.5) && (power < 2) ) {
      id.interp <- ( (xix > 0) & (xix < 0.9) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 1.5
        p.hi <- 2
        xix.lo <- 0
        xix.hi <- 0.9
        np <- 15
        nx <- 25
      }
    }
    
    ### Cases p>2   ###
    if ( (power > 2) && (power < 3) ) {
      id.interp <- ( (xix > 0) & (xix < 0.9) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 2
        p.hi <- 3
        xix.lo <- 0
        xix.hi <- 0.9
        np <- 15
        nx <- 25
      }
    }
    if ( (power >= 3) && (power < 4) ) {
      id.interp <- ( (xix > 0) & (xix < 0.9) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 3
        p.hi <- 4
        xix.lo <- 0
        xix.hi <- 0.9
        np <- 15
        nx <- 25
      }
    }
    if ( (power >= 4) && (power < 5) ) {
      id.interp <- ( (xix > 0) & (xix < 0.9) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 4
        p.hi <- 5
        xix.lo <- 0
        xix.hi <- 0.9
        np <- 15
        nx <- 25
      }
    }
    if ( (power >= 5) && (power < 7) ) {
      id.interp <- ( (xix > 0) & (xix < 0.5) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 5
        p.hi <- 7
        xix.lo <- 0
        xix.hi <- 0.5
        np <- 15
        nx <- 25
      }
    }
    if ( (power >= 7) && (power <= 10) ) {
      id.interp <- ( (xix > 0) & (xix < 0.3) )
      id.series <- (!(id.interp | id.type0))
      if ( any(id.interp)) {
        grid <- stored_grids(power)
        p.lo <- 7
        p.hi <- 10
        xix.lo <- 0
        xix.hi <- 0.3
        np <- 15
        nx <- 25
      }
    }
    
    if ( power > 10) {
      id.series <- (y != 0)
      id.interp <- (!(id.series | id.type0))
    }
    
    if (any(id.series)) {
      density[id.series] <- dtweedie_series(y = y[id.series],
                                            mu = mu[id.series], 
                                            phi = phi[id.series],
                                            power = power)
    }
    
    if (any(id.interp)) {
      dim( grid ) <- c( nx + 1, np + 1 )
      rho <- dtweedie_interp(grid, 
                             np = np, 
                             nx = nx,
                             xix.lo = xix.lo, 
                             xix.hi = xix.hi,
                             p.lo = p.lo, 
                             p.hi = p.hi,
                             power = power, 
                             xix = xix[id.interp] )
      dev <- tweedie_dev(power = power, 
                         mu = mu[id.interp], 
                         y = y[id.interp])
      front <- rho / (y[id.interp] * sqrt(2 * pi * xi[id.interp]))
      density[id.interp] <- front * exp(-1/(2 * phi[id.interp]) * dev)
    }
  }
  
  # Restoration of whole vector, and sanity fixes
  density2[ !special_y_cases ] <- density
  density <- density2
  if (any(density < 0 ) )  density[ density < 0 ] <- rep(0, sum(density < 0) )
  density <- as.vector(density)

  if (details) {
    return(list(density = density,
                regions = regions) )
  } else {
    return(density)
  }
  
}
