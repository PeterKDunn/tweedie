#############################################################################
ptweedie <- function(q, xi = NULL, mu, phi, power = NULL) {
  # Evaluates the cdf for Tweedie distributions, for given values of:
  #   q (possibly a vector)
  #   mu, the mean 
  #   phi the dispersion parameter
  #   power,  the Tweedie index parameter
  
  # Peter Dunn
  # Created: 01 May 2001
  # Last edit: 16 Sep 2025

    
  # SORT OUT THE NOTATION (i.e., xi VS power)
  out <- sort_Notation(xi = xi, 
                       power = power)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long ### MAY NOT BE NEEDED!!!
  
  # CHECK THE INPUTS ARE OK
  y <- q
  y.negative <- (y < 0)
  y.original <- y
  y <- y[ !y.negative ]
  
  out <- check_Inputs(y, mu, phi, power)
  mu <- out$mu
  phi <- out$phi
  mu.original <- out$mu.original
  phi.original <- out$phi.original
  
  
  # Set up arrays for CDF
  cdf.positives <- array( dim = length(y) )
  cdf <- y.original
  
  
  ### EVALUATE SPECIAL CASES: where p = 0, 1 or 2
  ### COULD WE ADD p = 3, based on statmod::pinvgauss
  out <- special_Cases_CDF(y, mu, phi, power)
  f <- out$f
  
  # Now, for p > 2 the only option is the inversion
  if ( power> 2 ) {
    
    # When y/q is very small, the function should return 0.
    # Sometimes it fails (when *very* small...).
    # This adjustment is made in ptweedie.inversion()
    f <- ptweedie.inversion(power = power, 
                            mu = mu, 
                            q = y, 
                            phi = phi)
  }
  
  # For 1<p<2, the two options are the series or inversion.
  # We avoid the series when p is near one, otherwise it is fine.
  
  # A few caveats.  Gustaov noted this case:
  # ptweedie(q=7.709933e-308, mu=1.017691e+01, phi=4.550000e+00, power=1.980000e+00)
  # which fails for the inversion, but seems to go fine for the series.
  # So the criterion should be a bit more detailed that just p<1.7...
  # But what?
  # Shallow second derivative of integrand is OK in principle... but hard to ascertain.
  
  # In a bug report by Gustavo Lacerda (April 2017), 
  # it seems that the original code here was a bit too 
  # harsh on the series, and a bit forgiving on the inversion.
  # Changed April 2017 to use the series more often.
  
  # ### OLD CODE:
  # if ( (power>1) & (power<2) ) {
  #     if ( power <1.7 ) {
  #        f <- ptweedie.series(power=power, q=y, mu=mu, phi=phi )
  #     } else{
  #        f <- ptweedie.inversion( power=power, q=y, mu=mu, phi=phi)
  #     }
  # }
  
  ### REVISED CODE:
  ### The choice of  1.999 is arbitrary.  Probably needs serious attention to decide properly
  ### Changed early 2017; thanks to Gustavo Lacerda
  if ( (power > 1) & (power < 2) ) {
    
    if ( power < 1.999) { 
      #### XXXXXXXXXXXXXXXXXXXXXXXXX This is arbitrary, and needs a closer look
      f <- ptweedie.series(power = power, 
                           q = y, 
                           mu = mu, 
                           phi = phi )
    } else{
      f <- ptweedie.inversion( power = power, 
                               q = y, 
                               mu = mu, 
                               phi = phi)
    }
  }
  
  # Sanity fixes
  cdf[ !y.negative ] <- f
  cdf[  y.negative ] <- 0
  cdf[ is.infinite( cdf ) ] <- 1
  cdf[ cdf > 1 ] <- 1
  
  return(cdf)
}

#############################################################################
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


#############################################################################
ptweedie.inversion <- function(q, mu, phi,  power, exact = FALSE ){ 
  # Evaluates the cdf for Tweedie distributions, using Fourier inversion, 
  # for given values of:
  #   q (possibly a vector)
  #   mu, the mean 
  #   phi the dispersion parameter
  #   power,  the Tweedie index parameter
  
  # Peter K Dunn 
  # Created: 08 Feb 2000 
  # Last edit: 16 Sep 2025

  y <- q
  cdf <- array( dim = length(y) )
  
  # Error checks
  if ( power < 1) stop("power must be greater than 1.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y < 0) ) stop("y must be a non-negative vector.")
  if ( any(mu <= 0) ) stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  #its <- y
  for (i in (1:length(y))) {
    
    # This has been added to avoid an issue with *very* small values of y
    # causing the FORTRAN to die when p>2; 
    # reported by Johann Cuenin 09 July 2013 (earlier, but that was the easiest email I could find about it :->)
    if ( ( power > 2 ) & (y[i] < 1.0e-300) ) {
      ### THE  e-300  IS ARBITRARY!!!!  ###
      # Keep an eye on it; perhaps it needs changing
      
      # That is, y is very small: Use the limit as y->0 as the answer as FORTRAN has difficulty converging
      # I have kept the call to the FORTRAN (and for p>2, of course, the limiting value is 0).
      #  I could have done this differently
      # by redefining very small y as y=0.... but this is better methinks
      cdf[i] <- 0
    } else {
      
      #   cat(power,phi[i], y[i], mu[i], as.integer(exact),"\n")
      tmp <- .Fortran( "twcdf",
                       as.double(power),
                       as.double(phi[i]),
                       as.double(y[i]),
                       as.double(mu[i]),
                       as.integer( exact ),
                       as.double(0), # funvalue
                       as.integer(0), # exitstatus
                       as.double(0), # relerr
                       as.integer(0)) # its
      cdf[i] <- tmp[[6]]
    }
  }
  
  cdf
  
}


