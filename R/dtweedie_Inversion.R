dtweedie.inversion <- function(y, power, mu, phi, method=3, verbose = FALSE, details = FALSE){ 
  # 
  # Peter K Dunn 
  # 06 Aug 2002
  #
  
  # Error checks
  if ( power < 1) stop("power must be greater than 1.")
  if ( any(phi <= 0)) stop("phi must be positive.")
  #if ( power>2) if ( any(y <= 0) ) stop("y must be a positive vector.")
  #if ( (power>1 & (power<2) ) if ( any(y < 0) ) stop("y must be a non-negative vector.")
  if ( any(mu <= 0) ) stop("mu must be positive.")
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
  
  N <- length(y)
  density <- y
  its <- y

  if ( is.null(method)){
    method <- array( dim = N)
  } else {
    method <- array( method, 
                     dim = N)
  }
  
  for (i in (1:N)) {
    
    # There are three approaches, each a product of a simple bit
    # and a complicated bit computed in FORTRAN
    # We choose Method 3 if no other is requested.
    #
    # The methods are documented in Dunn and Smyth (2008):
    # - Method 1: Evaluate a: compute a(y, phi) = f(y; 1, phi)
    # - Method 2: Rescale the mean to 1
    # - Method 3: Rescale y to 1 and evaluate b.
    
    if ( y[i] <= 0 ) {
      if ( (power > 1) & (power < 2) ) {
        if ( y[i] == 0 ) {
          density[i] <- exp( -mu[i] ^ (2 - power) / ( phi[i] * (2 - power) ) )
        } else {
          density[i] <- 0
        }
      } else {
        density[i] <- 0
      }
    } else {
      # Here, y > 0
      m2 <- 1 / mu[i]
      
      theta <- ( mu[i] ^ (1 - power) - 1 ) / ( 1 - power )
      if ( ( abs(power - 2 ) ) < 1.0e-07 ){
        kappa <- log(mu[i]) + (2 - power) * ( log(mu[i]) ^ 2 ) / 2
      } else {
        kappa <- ( mu[i] ^ (2 - power) - 1 ) / (2 - power)
      }
      m1 <- exp( (y[i]*theta - kappa ) / phi[i] )
      
      dev <- tweedie.dev(y = y[i], 
                         mu = mu[i],
                         power = power )
      m3 <- exp( -dev/(2 * phi[i]) ) / y[i]
      
      min.method <- min( m1, m2, m3 )
      
      # Now if no method requested, find the notional "optimal"
      if ( is.null(method[i]) ) {
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
        use.method <- method[i]
      }
      
      # Now use the method
      # NOTE: FOR ALL  METHODS, WE HAVE mu=1
      
      if ( use.method == 2 ) {
        tmp <- .C("twcdf",
                  N          = as.integer(N),
                  power      = as.double(power),
                  phi        = as.double(phi[i] / (mu[i] ^ (2 - power)) ), # phi
                  y          = as.double(y[i] / mu[i]),            # y
                  mu         = as.double(1),                       # mu
                  verbose    = as.integer(verbose),                # verbosity
                  pdf        = as.integer(1),                      # 1: TRUE. This is the PDF
                  # THE OUTPUTS:
                  funvalue   = as.double(0),                       # funvalue
                  exitstatus = as.integer(0),                      # exitstatus
                  relerr     = as.double(0),                       # relerr
                  its        = as.integer(0),                      # its
                  PACKAGE    = "tweedie")
        den <- tmp$funvalue
        density[i] <- den * m2
        
      } else {
        if ( use.method == 1 ) {
          tmp <- .C("twcdf",
                    N          = as.integer(N),
                    power      = as.double(power),
                    phi        = as.double(phi[i]),      # phi
                    y          = as.double(y[i]),        # y
                    mu         = as.double(1),           # mu = 1 for PDF
                    verbose    = as.integer( verbose ),  # verbose as an integer
                    pdf        = as.integer(1),          # 1: TRUE. This is the PDF
                    # THE OUTPUTS:
                    funvalue   = as.double(0),           # funvalue
                    exitstatus = as.integer(0),          # exitstatus
                    relerr     = as.double(0),           # relerr
                    its        = as.integer(0),          # its
                    PACKAGE    = "tweedie")
          
          den <- tmp[[7]]
          density[i] <- den * m1
          
        } else { # use.method == 3
          tmp <- .C("twcdf",
                    N          = as.integer(N),
                    power      = as.double(power),
                    phi        = as.double(phi[i] / (y[i] ^ (2 - power))), # phi
                    y          = as.double(1),          # y
                    mu         = as.double(1),          # mu
                    verbose    = as.integer( verbose ), # verbose as an integer
                    pdf        = as.integer(1),         # 1: TRUE. This is the PDF
                    # THE OUTPUTS:
                    funvalue   = as.double(0),          # funvalue
                    exitstatus = as.integer(0),         # exitstatus
                    relerr     = as.double(0),          # relerr
                    its        = as.integer(0),         # its
                    PACKAGE    = "tweedie")
          
          den <- tmp$funvalue
          density[i] <- den * m3
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

