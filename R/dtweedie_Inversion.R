dtweedie.inversion <- function(y, power, mu, phi, method = 3, verbose = FALSE, details = FALSE){ 
  # Evaluates the pdf for Tweedie distributions, using Fourier inversion, in FORTRAN:
  #
  #   y           : the values at which to compute the density (possibly a vector)
  #   power       : the Tweedie index parameter
  #   mu          : the mean (possibly a vector)
  #   phi         : the dispersion parameter (possibly a vector)
  #   method      : the method of evaluation used (as explained in Dunn & Smyth, 2008)
  #   verbose     : the verbosity of the output
  #   details     : whether to return details of the algorithm
  
  # Peter K Dunn 
  # Created: 08 Feb 2000 
  # Last edit: 28 Oct 2025
  
  
  ### BEGIN CHECKS
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) {
      stop("mu must be scalar, or the same length as y.")
    }
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
  save.method <- method
  if ( !is.null(method)){
    if ( length(method) > 1 ) {
      method <- save.method <- method[1]
    }
    if ( !(method %in% c(1, 2, 3)) ) stop("method must be 1, 2 or 3 (or left empty).")
  }
  ### END CHECKS
  

  ### BEGIN SET UP  
  N <- length(y)
  density <- y
  its <- y
  pSmall  <- ifelse( (power > 1) & (power < 2), TRUE, FALSE )

  if ( is.null(method)){
    method <- array( dim = N)
  } else {
    method <- array( method, 
                     dim = N)
  }
  ### END SET UP
  
  
  
  # There are three approaches ('method'), each a product of a simple bit
	# and a complicated bit computed in FORTRAN
  # We choose Method 3 if no other is requested.
  #
  # The methods are documented in Dunn and Smyth (2008):
  # - Method 1: Evaluate a(): compute a(y, phi) = f(y; 1, phi)
  # - Method 2: Rescale the mean to 1
  # - Method 3: Rescale y to 1 and evaluate b().
  
  
  if ( length(y0) > 0){
    if (pSmall) {
      lambda <- mu[y0]^(2 - power) / ( phi[y0] * (2 - power) )
      cdf[ y0 ] <- exp( -lambda )
    } else {
      cdf[ y0 ] <- rep(0,
                       length(y0) )
    }
  }
  NFortran <- N - length(y0)

 
 
 	### BEGIN ESTABLISH METHOD
	m2 <- 1 / mu[yFortran]
    
  theta <- ( mu[yFortran] ^ (1 - power) - 1 ) / ( 1 - power )
  if ( ( abs(power - 2 ) ) < 1.0e-07 ){
    kappa <- log(mu[yFortran]) + (2 - power) * ( log(mu[yFortran]) ^ 2 ) / 2
  } else {
    kappa <- ( mu[yFortran] ^ (2 - power) - 1 ) / (2 - power)
  }
  m1 <- exp( (y[yFortran] * theta - kappa ) / phi[yFortran] )
    
  dev <- tweedie.dev(y = y[yFortran], 
                     mu = mu[yFortran],
                     power = power )
  m3 <- exp( -dev/(2 * phi[yFortran]) ) / y[yFortran]
    
  min.method <- min( m1, m2, m3 )
    
  # Now if no method requested, find the notional "optimal"
  if ( is.null(method[yFortran]) ) {
    if ( min.method == m1 ){
      use.method <- 1
    } else {
      if ( min.method == m2 ) {
      	use.method <- 2
      } else {
      	use.method <- 3
      }
    }
  } else {
    use.method <- method[yFortran]
  }
  ## END ESTABLISH METHOD
    
    
    
  # Now use the method
  # NOTE: FOR ALL  METHODS, WE HAVE mu=1
	if (NFortran > 0 ) {
   
		if ( use.method == 2 ) {
			tmp <- .C("twcomputation",
						N          = as.integer(NFortran),
						power      = as.double(power),
						phi        = as.double(phi[yFortran] / (mu[yFortran] ^ (2 - power)) ), # phi
						y          = as.double(y[yFortran] / mu[yFortran]),            # y
						mu         = as.double(rep(1, NFortran)),                       # mu
						verbose    = as.integer(verbose),                # verbosity
						pdf        = as.integer(1),                      # 1: TRUE, as this is the PDF
						# THE OUTPUTS:
						funvalue   = as.double(rep(0, NFortran)),                       # funvalue
						exitstatus = as.integer(0),                      # exitstatus
						relerr     = as.double(0),                       # relerr
						its        = as.integer(rep(0, NFortran)),                      # its
						PACKAGE    = "tweedie")
			den <- tmp$funvalue
			density[yFortran] <- den * m2
			
			} else {
			if ( use.method == 1 ) {
				tmp <- .C("twcomputation",
						N          = as.integer(NFortran),
						power      = as.double(power),
						phi        = as.double(phi[yFortran]),      # phi
						y          = as.double(y[yFortran]),        # y
						mu         = as.double(rep(1, NFortran)),           # mu = 1 for PDF
						verbose    = as.integer( verbose ),  # verbose as an integer
						pdf        = as.integer(1),          # 1: TRUE, as this is the PDF
						# THE OUTPUTS:
						funvalue   = as.double(rep(0, NFortran)),           # funvalue
						exitstatus = as.integer(0),          # exitstatus
						relerr     = as.double(0),           # relerr
						its        = as.integer(rep(0, NFortran)),          # its
						PACKAGE    = "tweedie")
				
				den <- tmp[[7]]
				density[yFortran] <- den * m1
				
			} else { # use.method == 3
				tmp <- .C("twcomputation",
						N          = as.integer(NFortran),
						power      = as.double(power),
						phi        = as.double(phi[yFortran] / (y[yFortran] ^ (2 - power))), # phi
						y          = as.double(rep(1, NFortran)),          # y
						mu         = as.double(rep(1, NFortran)),          # mu
						verbose    = as.integer( verbose ), # verbose as an integer
						pdf        = as.integer(1),         # 1: TRUE, as this is the PDF
						# THE OUTPUTS:
						funvalue   = as.double(rep(0, NFortran)),          # funvalue
						exitstatus = as.integer(0),         # exitstatus
						relerr     = as.double(0),          # relerr
						its        = as.integer(rep(0, NFortran)),         # its
						PACKAGE    = "tweedie")
				
				den <- tmp$funvalue
				density[yFortran] <- den * m3
			}
    }
  }
  
  if (details) {
    return( list( density = den,
                  regions = tmp$its))
  } else {
    return(density)
  }
  
}

