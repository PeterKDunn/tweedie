#' Fourier Inversion Evaluation for the Tweedie Distribution Function
#'
#' @description
#' Evaluates the distribution function (\acronym{df}) for Tweedie distributions using Fourier inversion, 
#' for given values of the dependent variable \code{y}, 
#' the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be in the case of evaluation problems.
#'
#' @usage ptweedie_inversion(q, mu, phi, power, verbose = FALSE, details = FALSE)
#'
#' @param q vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter.
#' @param phi the dispersion parameter.
#' @param verbose logical; if \code{TRUE}, displays some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution and some information about the integration. The default is \code{FALSE}.
#' 
#' @return if \code{details = FALSE}, a numeric vector of the distribution function values; if \code{details = TRUE}, a list containing \code{CDF} (a vector of the values of the distribution function) and \code{regions} (a vector of the number of integration regions used).
#' 
#' @aliases ptweedie.inversion
#' @export
ptweedie_inversion <- function(q, mu, phi, power, verbose = FALSE, details = FALSE ){ 
  ### NOTE: No checking of inputs
  ### Assumes that all of y, mu, phi have the same length (they may be vectors) and are valid
  
  ### BEGIN SET UP
  N <- as.integer( length(q) )
  cdf <- as.double(rep(0, N))
  pSmall  <- ifelse( (power > 1) & (power < 2),
                     TRUE, 
                     FALSE )

  # Initialise
  exitstatus_scalar <- as.integer(0)
  relerr_scalar     <- as.double(0.0)
  its_scalar        <- as.integer(0)
  ### END SET UP


  tmp <- .C( "twcomputation",
             N           = as.integer(N),         # number of observations
             power       = as.double(power),      # p
             phi         = as.double(phi),        # phi
             y           = as.double(q),          # y
             mu          = as.double(mu),         # mu
             verbose     = as.integer(verbose),   # verbosity
             pdf         = as.integer(0),         # 0: FALSE, as this is the CDF not PDF
             # THE OUTPUTS:
             funvalue    = as.double(rep(0, N)),  # funvalue
             exitstatus  = as.integer(0),         # exitstatus
             relerr      = as.double(0),          # relerr
             its         = as.integer(rep(0, N)), # its
             PACKAGE     = "tweedie")
  cdf <- tmp$funvalue

  if (details) {
    return( list( cdf = cdf,
                  regions = tmp$its))
  } else {
    return(cdf)
  }
}

#' @export
ptweedie.inversion <- function(q, power, mu, phi, verbose, details){ 
  .Deprecated("ptweedie_inversion", package = "tweedie")
  ptweedie_inversion(q = q, 
                     power = power,
                     mu = mu, 
                     phi = phi, 
                     verbose = FALSE, 
                     details = FALSE)
}

