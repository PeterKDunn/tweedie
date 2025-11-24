#' Density for the Tweedie Family of Distributions
#'
#' Evaluates the AIC for a fitted Tweedie glm.
#' 
#' @description Evaluates the AIC for a fitted Tweedie glm.
#'   The Tweedie family of distributions belong to the class of
#'   exponential dispersion models (EDMs), famous for their role in generalized
#'   linear models. The Tweedie distributions are the EDMs with a variance
#'   of the form \eqn{\mbox{var}[Y] = \phi\mu^p}{var[Y] = phi*mu^p} where \eqn{p \ge 1}{p >= 1}.
#'   \bold{This function only evaluates for \eqn{p \ge 1}{p >= 1}.}
#'
#' @details
#' Evaluating the likelihood can be time consuming, so the function may take some time.
#'
#' @param glm.obj A fitted \code{glm} object, fitted using the \code{tweedie} family.
#' @param dispersion The dispersion parameter, usually extracted from \code{glm.obj}; however, occasionally  a specified value of the dispersion may be needed.
#' @param k The AIC penalty; \code{k = 2} (the default) produces the AIC.
#' @param verbose Logical. Whether to display details of the internal process. Defaults to \code{FALSE}.
#'
#' @return The value of the computed AIC.
#'
#' @seealso \code{\link{ptweedie}}, \code{\link{rtweedie}}, \code{\link{dtweedie}}
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
#' Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986). 
#' \emph{Akaike Information Criterion Statistics}. 
#' D. Reidel Publishing Company.
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
#' @aliases tweedie_AIC
#' @aliases AICtweedie
#' 
#' @export
tweedie_AIC <- function( glm.obj, dispersion = NULL, k = 2, verbose = TRUE){ 
  # New  dispersion  input for (e.g.) Poisson case, added October 2017
  
  wt <- glm.obj$prior.weights
  n <- length(glm.obj$residuals)
  edf <- glm.obj$rank  # As used in logLik.glm()
  
  mu <- stats::fitted( glm.obj )
  y  <- glm.obj$y
  p <- get("p", envir = environment(glm.obj$family$variance))
  
  if ( is.null(dispersion)) {  # New section
    if (p == 1 & verbose) message("*** Tweedie index power = 1: Consider using  dispersion=1  in call to  AICtweedie().\n")
    dev <- deviance(glm.obj)
    disp <- dev / sum(wt)  # In line with Gamma()$aic
    edf <- edf + 1  # ADD one as we are estimating phi too
  } else {
    disp <- dispersion
  }
  
  den <- dtweedie( y = y, 
                   mu = mu, 
                   phi = disp, 
                   power = p)
  AIC <- -2 * sum( log(den) * wt) 
  
  return( AIC + k * (edf) )
  
}




#' @export
AICtweedie <- function( glm.obj, dispersion = NULL, k = 2, verbose = TRUE){ 
  .Deprecated("tweedie_AIC", package = "tweedie")
  tweedie_AIC(glm.obj = glm.obj, 
              dispersion = dispersion,
              k = k,
              verbose = verbose)
}


