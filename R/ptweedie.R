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

  # SORT OUT THE NOTATION (i.e., xi VS power)
  if (verbose) cat("- Checking notation\n")
  out <- sort_Notation(xi = xi, power = power)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long ### MAY NOT BE NEEDED!!!
  
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_Inputs(q, mu, phi, power)
  mu <- out$mu
  phi <- out$phi
  f <- array(0,
             dim = length(q) )
  if (details) regions <- array(0, dim = length(q))
  

  # IDENTIFY SPECIAL CASES
  if (verbose) cat("- Checking for special cases\n")
  out <- special_Cases(q, mu, phi, power)
  f <- out$f
  special_p_Cases <- out$special_p_Cases
  if (verbose & special_p_Cases) cat("  - Special case for p used\n")
  special_y_Cases <- out$special_y_Cases
  if (verbose & any(special_y_Cases)) cat("  - Special cases for first input found\n")
  

  if ( !special_p_Cases ) {
    # NO special p cases
    if ( power > 2 ) {
      # For p > 2 the only option is the inversion
      if (verbose) cat("- With p > 2: use inversion")

      f_TMP <- ptweedie.inversion(power = power, 
                                  q = q[!special_y_Cases], 
                                  mu = mu[!special_y_Cases], 
                                  phi = phi[!special_y_Cases],
                                  verbose = verbose,
                                  details = details)
      if (details) {
        f[!special_y_Cases] <- f_TMP$cdf
        regions[!special_y_Cases] <- f_TMP$regions
      } else {
        f[!special_y_Cases] <- f_TMP
      }
    } else {
      # NOT a special-p case
      # There still may be some special-y cases
      
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
      if ( power < 1.999) { 
        #### XXXXXXXXXXXXXXXXXXXXXXXXX This is arbitrary, and needs a closer look
        if (verbose) cat("- With 1 < p < 2: use series")
        
        f_TMP <- ptweedie.series(power = power, 
                                 q     = q[!special_y_Cases], 
                                 mu    = mu[!special_y_Cases], 
                                 phi   = phi[!special_y_Cases] )
        if (details) {
          f[!special_y_Cases] <- f_TMP$cdf
          regions[!special_y_Cases] <- f_TMP$regions
        } else {
          f[!special_y_Cases] <- f_TMP
        }
      } else{
        if (verbose) cat("- With 1 < p < 2: use inversion")
        
        f_TMP <- ptweedie.inversion( power   = power, 
                                    q       = q[!special_y_Cases], 
                                    mu      = mu[!special_y_Cases], 
                                    phi     = phi[!special_y_Cases],
                                    verbose = verbose,
                                    details = details)
        if (details) {
          f[!special_y_Cases] <- f_TMP$cdf
          regions[!special_y_Cases] <- f_TMP$regions
        } else {
          f[!special_y_Cases] <- f_TMP
        }
      }
    }
  }
  
  
  # Sanity fixes
  f[ f < 0 ] <- 0
  f[ f > 1 ] < 1

  if (details) {
    return(list(f = f,
                regions = regions) )
  } else {
    return(f)
  }
}

