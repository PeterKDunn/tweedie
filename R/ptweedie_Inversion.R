ptweedie.inversion <- function(q, mu, phi,  power, verbose = FALSE, details = FALSE ){ 
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
  
  N <- as.integer( length(y) )
  cdf <- as.double(rep(0, N))
  
  small_Y_cases <- which( ( power > 2 ) & (y < 1.0e-300) )
  
  # Set CDF to 0 for the edge cases
  if (length(small_Y_cases) > 0) {
    cdf[small_Y_cases] <- 0
  }
  
  # Identify the non-edge cases that need Fortran calculation
  need_Fortran <- which( !(( power > 2 ) & (y < 1.0e-300)) )
  # This has been added to avoid an issue with *very* small values of y
  # causing the FORTRAN to die when p>2; 
  # reported by Johann Cuenin 09 July 2013 (earlier, but that was the easiest email I could find about it :->)
  ### THE  e-300  IS ARBITRARY!!!!  ###
  # Keep an eye on it; perhaps it needs changing
  # That is, y is very small: Use the limit as y->0 as the answer as FORTRAN has difficulty converging
  # I have kept the call to the FORTRAN (and for p>2, of course, the limiting value is 0).
  #  I could have done this differently
  # by redefining very small y as y=0.... but this is better methinks
  
  # If there's nothing to calculate, return
  if (length(need_Fortran) == 0) {
    return(cdf)
  }
  
  # Initialise
  exitstatus_scalar <- as.integer(0)
  relerr_scalar     <- as.double(0.0)
  its_scalar        <- as.integer(0)
  
  tmp <- .C( "twcomputation",
             N           = as.integer(N),         # number of observations
             power       = as.double(power),      # p
             phi         = as.double(phi),        # phi
             y           = as.double(y),          # y
             mu          = as.double(mu),         # mu
             verbose     = as.integer(verbose),   # verbosity
             pdf         = as.integer(0),         # 0: FALSE. This is the CDF not PDF
             # THE OUTPUTS:
             funvalue    = as.double(rep(0, N)),  # funvalue
             exitstatus  = as.integer(0),         # exitstatus
             relerr      = as.double(0),          # relerr
             its         = as.integer(0),         # its
             PACKAGE     = "tweedie")
  cdf <- tmp$funvalue
  
  if (details) {
    return( list( cdf = cdf,
                  regions = tmp$its))
  } else {
    return(cdf)
  }
}

