#' Tweedie Distributions: Plot Tweedie density functions
#'
#' Plot the Tweedie PDF.
#' @description This function produced a plot of the specified Tweedie distribution.
#'
#' @details If \eqn{1 < p < 2}{1 < power < 2}, the mass at \eqn{Y=0}{Y = 0} is automatically added.
#'
#' @param y The values for \eqn{y}{y} in the plot.
#' 
#' @param power The variance power.
#' 
#' @param mu The mean of the distribution.
#' 
#' @param phi The dispersion parameter.
#' 
#' @param xi A synonym for \code{power}.
#' 
#' @param type The type of plot, either \code{PDF} (the default) or \code{CDF}.
#'
#' @param add If \code{TRUE}, the plot is added to the current plot; if \code{FALSE} (the default) the plot is produced on a fresh plot.
#' 
#' @param ... Plotting parameters
#' 
#' @importFrom graphics lines rug par mtext abline axis  points plot
tweedie_Plot <- function(y, xi = NULL, mu, phi, type = "pdf", power = NULL, 
                         add =FALSE, ...) {
  
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
  
  if ( ( power < 0 ) | ( ( power > 0 ) & ( power < 1 ) ) ) {
    stop( paste("Plots cannot be produced for", index.par.long, "=", power, "\n") )
  }   
  
  is.pg <- ( power > 1 ) & ( power < 2 )
  
  if ( type == "pdf") {
    fy <- dtweedie( y = y, 
                    power = power, 
                    mu = mu,
                    phi = phi)
  } else {
    fy <- ptweedie( q = y, 
                    power = power, 
                    mu = mu,
                    phi = phi)
  }
  
  if ( !add ) {
    if ( is.pg ) {
      plot( range(fy) ~ range( y ), 
            type = "n", ...)
      if ( any( y == 0 ) ) { # The exact zero
        points( fy[y == 0] ~ y[y == 0], 
                pch = 19, ... )
      }
      if ( any( y > 0 ) ) { # The exact zero
        lines( fy[y > 0]   ~ y[y > 0], 
               pch = 19, ... )
      }
    } else {  # Not a Poison-gamma dist
      plot( range(fy) ~ range( y ), 
            type = "n", ...)
      lines( fy ~ y, 
             pch = 19, ... )
    }
  } else {# Add; no new plot
    if ( is.pg ) {
      if ( any( y == 0 ) ) { # The exact zero
        points( fy[y == 0] ~ y[y == 0], 
                pch = 19, ... )
      }
      if ( any( y > 0 ) ) { # The exact zero
        lines( fy[y > 0]   ~ y[y > 0], 
               pch = 19, ... )
      }
    } else {  # Not a Poison-gamma dist
      lines( fy ~ y, 
             pch = 19, ... )
    }
    
  }
  return(invisible(list(y = fy, 
                        x = y) ))
  
}