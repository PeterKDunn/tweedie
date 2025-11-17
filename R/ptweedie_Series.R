#' Tweedie Distribution Function: Series Evaluation for the Distribution Function
#'
#' Internal function to evaluate the Tweedie distribution function using the infinite series expansion.
#' \bold{Not intended for general users.}
#'
#' @param q vector of quantiles.
#' @param power The power parameter \eqn{p}{power}.
#' @param mu The mean parameter.
#' @param phi The dispersion parameter.
#' @param verbose Display some internal computation details.
#' @param details Return the DF and the number of series terms used.
#' @return A numeric vector of densities.
#' 
#' @importFrom stats dpois 
#'
#' @keywords internal

ptweedie.series <- function(q, power, mu, phi, verbose = FALSE, details = FALSE) {

  # Evaluates the cdf for Tweedie distributions, using the infinite series 
  # expansion (i.e., just for 1 < p < 2), for given values of:
  #   q (possibly a vector)
  #   mu, the mean 
  #   phi the dispersion parameter
  #   power,  the Tweedie index parameter
  
  ### NOTE: No checking of inputs
  ### Assumes that all of y, mu, phi have the same length (they may be vectors) and are valid
  
  ### NOTE: Same number of iterations used for each value of q, if vector inputs are given

  # SET UP
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  drop <- 39

  # FIND THE LIMITS ON N, the summation index
  # The *lower* limit on N
  lambda <- max(lambda )
  logfmax <-  -log(lambda)/2
  estlogf <- logfmax
  N <- max( lambda )
  
  while ( ( estlogf > (logfmax - drop) ) & ( N > 1 ) ) {
    N <- max(1, N - 2)
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  lo.N <- max(1, floor(N) )
  
  
  # The *upper* limit on N
  lambda <- min( lambda )
  logfmax <-  -log(lambda) / 2
  estlogf <- logfmax
  N <- max( lambda )
  
  while ( estlogf > (logfmax - drop) ) {
    N <- N + 1
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  hi.N <- max( ceiling(N) )
  if (verbose) cat("Summing over", lo.N, "to", hi.N, "\n")
  
  # EVALUATE between limits of N
  cdf <- array( dim = length(q), 0 )
  
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  
  
  for (N in (lo.N : hi.N)) {
    # Poisson density
    pois.den <- dpois( N, lambda)
    
    # Incomplete gamma
    incgamma.den <- stats::pchisq(2 * q / tau, 
                           -2 * alpha * N )
    
    # What we want
    cdf <- cdf + pois.den * incgamma.den
    
  }
  
  cdf <- cdf + exp( -lambda )
  its <- hi.N - lo.N + 1
  
  if (details) {
    return( list( cdf = cdf,
                  iterations = its) )
  } else {
    return(cdf)
  }
  
}

