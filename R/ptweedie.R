#############################################################################
ptweedie <- function(q, xi = NULL, mu, phi, power = NULL, verbose = FALSE, details = FALSE) {
  # Evaluates the cdf for Tweedie distributions, for given values of:
  #   q (possibly a vector)
  #   mu, the mean 
  #   phi the dispersion parameter
  #   power,  the Tweedie index parameter
  #   verbose: whether to display what is happening
  #   details:  whether to returns reports of relerr, regions needed, etc.
  
  # Peter Dunn
  # Created: 01 May 2001
  # Last edit: 16 Sep 2025

cat(">> ptweedie: GOT THIS FAR\n")    
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
                            phi = phi,
                            verbose = verbose,
                            details = details)
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
      f <- ptweedie.inversion( power   = power, 
                               q       = y, 
                               mu      = mu, 
                               phi     = phi,
                               verbose = verbose,
                               details  = details)
    }
  }
  
  # Sanity fixes
  cdf[ !y.negative ] <- f
  cdf[  y.negative ] <- 0
  cdf[ is.infinite( cdf ) ] <- 1
  cdf[ cdf > 1 ] <- 1
  
  return(cdf)
}

