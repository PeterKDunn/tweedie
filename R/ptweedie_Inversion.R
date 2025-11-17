#' Tweedie Distribution: Fourier Inversion Evaluation for the Distribution Function
#'
#' Internal function to evaluate the Tweedie distribution function using Fourier inversion.
#' \bold{Not intended for general users.}
#'
#' @param q vector of quantiles.
#' @param power The power parameter \eqn{p}{power}.
#' @param mu The mean parameter.
#' @param phi The dispersion parameter.
#' @param verbose Display some internal computation details.
#' @param details Return the DF and the number of integral regions used.
#' @return A numeric vector of the distribution function values
#' @keywords internal
#' @export
ptweedie.inversion <- function(q, mu, phi, power, verbose = FALSE, details = FALSE ){ 
  # Evaluates the pdf for Tweedie distributions, using Fourier inversion, in FORTRAN:
  #
  #   q           : the values at which to compute the DF (possibly a vector)
  #   power       : the Tweedie index parameter
  #   mu          : the mean (possibly a vector)
  #   phi         : the dispersion parameter (possibly a vector)
  #   verbose     : the verbosity of the output
  #   details     : whether to return details of the algorithm
  

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

