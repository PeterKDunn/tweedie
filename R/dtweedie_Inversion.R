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
  
  ### NOTE: No checking of inputs
  
  ### BEGIN SET UP
  N <- as.integer( length(y) )
  density <- as.double(rep(0, N))
  pSmall  <- ifelse( (power > 1) & (power < 2), TRUE, FALSE )

  # Initialise
  exitstatus_scalar <- as.integer(0)
  relerr_scalar     <- as.double(0.0)
  its_scalar        <- as.integer(0)
  ### END SET UP
  
  
  # Establish which method to use
  if ( is.null(method)){
    method <- array( dim = N)
  } else {
    method <- array( method, 
                     dim = N)
  }
  # There are three approaches ('method'), each a product of a simple bit
	# and a complicated bit computed in FORTRAN
  #
  # The methods are documented in Dunn and Smyth (2008):
  # - Method 1: Evaluate a(): compute a(y, phi) = f(y; 1, phi)
  # - Method 2: Rescale the mean to 1
  # - Method 3: Rescale y to 1 and evaluate b().
  # If no method is explicitly requested, find the notional "optimal" method for each i.
  

 	### BEGIN ESTABLISH METHOD
	m2 <- 1 / mu
    
  theta <- ( mu ^ (1 - power) - 1 ) / ( 1 - power )
  if ( ( abs(power - 2 ) ) < 1.0e-07 ){
    kappa <- log(mu) + (2 - power) * ( log(mu) ^ 2 ) / 2
  } else {
    kappa <- ( mu ^ (2 - power) - 1 ) / (2 - power)
  }
  m1 <- exp( (y * theta - kappa ) / phi )
    
  dev <- tweedie.dev(y = y, 
                     mu = mu,
                     power = power )
  m3 <- exp( -dev/(2 * phi) ) / y
    
  min.method <- min( m1, m2, m3 )
    
  # If no method is explicitly requested, find the notional "optimal" method for each i.
  use.method <- array( 3L, length = length(y) )
  if ( is.null(method) ) {
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
    use.method <- method
  }
  ## END ESTABLISH METHOD
    
  
  ### BEGIN COMPUTE

  
  ### END COMPUTE
  
  
  
  
  
  if (FALSE) {
    
    
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
}
  
  
  
  
  
  
  if (details) {
    return( list( density = den,
                  regions = tmp$its))
  } else {
    return(density)
  }
  
}

