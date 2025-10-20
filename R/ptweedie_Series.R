ptweedie.series <- function(q, power, mu, phi) {
  # Evaluates the cdf for Tweedie distributions, using the infinite series 
  # expansion (i.e., just for 1 < p < 2), for given values of:
  #   q (possibly a vector)
  #   mu, the mean 
  #   phi the dispersion parameter
  #   power,  the Tweedie index parameter
  
  # Peter K Dunn 
  # Created: 14 March 2000 
  # Last edit: 16 Sep 2025
  
  y <- q
  
  # Error checks
  if ( power < 1) stop("power must be between 1 and 2.\n")
  if ( power > 2) stop("power must be between 1 and 2.\n")
  if ( any(phi <= 0)) stop("phi must be positive.\n")
  if ( any(y < 0) ) stop("y must be a non-negative vector.\n")
  if ( any(mu <= 0) ) stop("mu must be positive.\n")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.\n")
  } else {
    mu <- array( dim = length(y), mu )
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.\n")
  } else {
    phi <- array( dim = length(y), phi )
  }
  
  
  # Set up
  p <- power
  lambda <- mu ^ (2 - p) / ( phi * (2 - p) )
  tau <- phi * (p - 1) * mu ^ ( p - 1 )
  alpha <- (2 - p) / (1 - p)
  
  
  # Now find the limits on N, the summation index
  # First the *lower* limit on N
  
  drop <- 39
  lambda <- max(lambda )
  logfmax <-  -log(lambda)/2
  estlogf <- logfmax
  N <- max( lambda )
  
  while ( ( estlogf > (logfmax - drop) ) & ( N > 1 ) ) {
    N <- max(1, N - 2)
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  
  lo.N <- max(1, floor(N) )
  
  
  # Now for the *upper* limit on N
  
  lambda <- mu ^ (2 - p) / ( phi * (2 - p) )
  lambda <- min( lambda )
  logfmax <-  -log(lambda) / 2
  estlogf <- logfmax
  
  N <- max( lambda )
  
  while ( estlogf > (logfmax - drop) ) {
    N <- N + 1
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  
  hi.N <- max( ceiling(N) )
  
  
  # Now evaluate between limits of N
  
  cdf <- array( dim = length(y), 0 )
  
  lambda <- mu ^ (2 - p) / ( phi * (2 - p) )
  tau <- phi * (p - 1) * mu ^ ( p - 1 )
  alpha <- (2 - p) / (1 - p)
  
  
  for (N in (lo.N : hi.N)) {
    
    # Poisson density
    pois.den <- dpois( N, lambda)
    
    # Incomplete gamma
    incgamma.den <- pchisq(2 * y / tau, 
                           -2 * alpha * N )
    
    # What we want
    cdf <- cdf + pois.den * incgamma.den
    
  }
  
  cdf <- cdf + exp( -lambda )
  its <- hi.N - lo.N + 1
  
  cdf
  
}

