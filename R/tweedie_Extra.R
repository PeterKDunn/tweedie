ptweedie_Inversion_Report <- function(q, mu, phi,  power, exact=FALSE ){ 
  #
  # ptweedie_INVERSION_Report
  #
  # ptweedie.inversion( y, power, mu, phi, exact )
  #
  # Evaluating Tweedie cdf using cgf inversion
  #
  # Arguments:
  #   y:  The values of a the random variable.  y  can be 
  #       a vector, and all elements must be non-negative.
  #
  #   mu: The (positive) mean of the distribution.
  #
  #   phi:    The (positive) dispersion parameter for the 
  #       distribution, such that the variance of the 
  #       observations is  var(Y) = phi * mu^power
  #
  #   power:  The index for the Tweedie distribution, such that 
  #       the variance of the observation is 
  #       var(Y) = phi * mu^power .  
  #
  #   exact:  If exact is TRUE, the exact zero acceleration algorithm 
  #       is used, otherwise the approx zero algorithm is used
  #
  # Description
  #   This function evaluates the Tweedie family of distributions
  #   distribution function by inverting the cgf.
  #
  # Value
  #   A list is returned with the following components:
  #   density:    The value of the density 
  #   
  # See Also:
  #   dtweedie.series.bigp (using series evaluation for evaluation)
  #   dtweedie (density evaluation using fast methods based on the
  #       series evaluation and cgf inversion)
  #
  # Examples
  #   dtweedie.inversion.report(mu=1, phi=1, y=data, power=3.3)
  
  
  # 
  # Peter K Dunn 
  # 08 Feb 2000 
  #
  
  y <- q
  # Error checks
  if ( power < 1) stop("power must be greater than 1...")
  if ( any(phi <= 0)) stop("phi must be positive.")
  if ( any(y < 0) ) stop("y must be a non-negative positive vector.")
  if ( any(mu <= 0) ) stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  }
  else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi)>1) {
    if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
  }
  else {
    phi <- array( dim=length(y), phi )
    # A vector of all phi's
  }
  
  
  y.len <- length(y)
  density <- y
  its <- y
  
  check_Inputs(y, mu, phi, power)
  
  
  for (i in (1:y.len)) {
    
    # This has been added to avoid an issue with *very* small values of y
    # causing the FORTRAN to die when p>2; 
    # reported by Johann Cuenin 09 July 2013 (earlier, but that was the easiest email I could find about it :->)
    if ( ( power > 2 ) & (y[i] < 1.0e-300) ) {
      # That is, y is very small: Use the limit as y->0 as the answer as FORTRAN has difficulty converging
      # I have kept the call to the FORTRAN (and for p>2, of course, the limiting value is 0).
      #  I could have done this differently
      # by redefining very small y as y=0.... but this is better methinks
      density[i] <- 0
      #its[i] <- tmp[[9]]
    } 
    else{
      tmp <- .Fortran( "twcdf", 
                       as.double(power),
                       as.double(phi[i]),
                       as.double(y[i]),
                       as.double(mu[i]),
                       as.double(0), # funvalue
                       as.integer(0), # exitstatus
                       as.double(0), # relerr
                       as.integer(0), # its
                       PACKAGE="tweedie")
      
      density[i] <- tmp[[6]]
      its[i] <- tmp[[9]]
    }
  }
  
  list(density = density, 
       its = its)
  
}



dtweedie_Series_Report <- function(y, power, mu, phi){ 
  #
  # dtweedie_SERIES_Report
  #
  # Evaluating Tweedie densities using series evaluation
  #
  # dtweedie.series( y, power, mu, phi)
  #
  # Arguments:
  #	power:  The index for the Tweedie distribution, such that 
  #		the variance of the observation is 
  #		var(Y) = phi * mu^power .
  #		Here, we require power>1
  #
  #	y:	The values of the random variable.  y  can be a 
  #		vector,	and all elements must be strictly positive.
  #
  #	mu:	The (positive) mean of the distribution
  #
  #	phi:	The (positive) dispersion parameter for the 
  #		distribution, such that the variance of the 
  #		observations is  var(Y) = phi * mu^power .
  #
  # Description
  #	This function evaluates the Tweedie family of distributions
  #	by evaluating an infinite series.
  #
  # Value
  #	The value of the density at the given parameters
  #	
  # See Also:
  #	dtweedie.series.bigp (for use when power > 2 )
  #	dtweedie.series.smallp (for use when 1 < power < 2 )
  #	dtweedie (density evaluation using fast methods based on the
  #		series evaluation and cgf inversion)
  #
  # Examples
  # 	dtweedie.series(mu=1, phi=1, y=data, power=1.3)
  
  
  #
  # Peter K Dunn
  # 09 Jan 2002
  #
  
  #
  # Error traps
  #
  
  if ( power<1) stop("power must be between 1 and 2.")
  if ( any(phi<=0) ) stop("phi must be positive.")
  if ( any(y<0) ) stop("y must be a non-negative vector.")
  if ( any(mu<=0) ) stop("mu must be positive.")
  if ( length(mu)>1) {
    if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
  }
  else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi)>1) {
    if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
  }
  else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  lo <- hi <- array( dim=length(y), NA )
  y0 <- (y == 0 )
  yp <- ( y!=0 )
  density <- y
  
  if ( (power == 2) | (power==1) ) { # Special cases
    if ( power == 2 ){
      density <- dgamma( y, shape=1/phi, rate=1/(phi * mu ) )
    }
    if ( power == 1 ){
      density <- dpois( y, lambda=mu )
    }
  }
  else{
    
    if ( any(y==0) ) {
      if ( power>2 ) {
        density[y0] <- 0*y[y0]
      }
      if ( ( power>1) && (power<2) ) {
        lambda <- mu[y0]^(2-power) / ( phi[y0] * (2-power) )
        density[y0] <- exp( -lambda )
      }
    }
    
    if ( any( y!=0 ) ) { 
      if ( power > 2 ) {
        out <- dtweedie.series.bigp(
          power=power,mu=mu[yp], y=y[yp], phi=phi[yp])
        density[yp] <- out$density
        lo[yp] <- out$lo
        hi[yp] <- out$hi
      }
      
      if ( ( power > 1 ) && ( power < 2 ) ) {
        out <- dtweedie.series.smallp(
          power=power,mu=mu[yp], y=y[yp], phi=phi[yp])
        density[yp] <- out$density
        lo[yp] <- out$lo
        hi[yp] <- out$hi
      }
    }
  }
  
  list(density = density, 
       lo = lo, 
       hi = hi)
  
}



#ptweedie.inversion.report <- function(q, mu, phi,  power, exact=FALSE ){ 
#
# ptweedie.INVERSION.report
#
# ptweedie.inversion( y, power, mu, phi, exact )
#
# Evaluating Tweedie cdf using cgf inversion
#
# Arguments:
#   y:  The values of a the random variable.  y  can be 
#       a vector, and all elements must be non-negative.
#
#   mu: The (positive) mean of the distribution.
#
#   phi:    The (positive) dispersion parameter for the 
#       distribution, such that the variance of the 
#       observations is  var(Y) = phi * mu^power
#
#   power:  The index for the Tweedie distribution, such that 
#       the variance of the observation is 
#       var(Y) = phi * mu^power .  
#
#   exact:  If exact is TRUE, the exact zero acceleration algorithm 
#       is used, otherwise the approx zero algorithm is used
#
# Description
#   This function evaluates the Tweedie family of distributions
#   distribution function by inverting the cgf.
#
# Value
#   A list is returned with the following components:
#   density:    The value of the density 
#   
# See Also:
#   dtweedie.series.bigp (using series evaluation for evaluation)
#   dtweedie (density evaluation using fast methods based on the
#       series evaluation and cgf inversion)
#
## Examples
##   dtweedie.inversion(mu=1, phi=1, y=data, power=3.3)
# 
# 
## 
## Peter K Dunn 
## 08 Feb 2000 
##
#
#y <- q
## Error checks
#if ( power<1) stop("power must be greater than 1.")
#if ( any(phi)<= 0) stop("phi must be positive.")
#if ( any(y<0) ) stop("y must be a non-negative positive vector.")
#if ( any(mu)<=0 ) stop("mu must be positive.")
#if ( length(mu)>1) {
#   if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
#}
#else {
#   mu <- array( dim=length(y), mu )
#   # A vector of all mu's
#}
#if ( length(phi)>1) {
#   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
#}
#else {
#   phi <- array( dim=length(y), phi )
#   # A vector of all phi's
#}
#   
#
#y.len <- length(y)
#density <- y
#its <- y
#
#for (i in (1:y.len)) {
#
#   tmp <- .Fortran( "cdf",
#          as.double(power),
#          as.double(phi[i]),
#          as.double(y[i]),
#          as.double(mu[i]),
#          as.integer( exact ),
#          as.double(0), # funvalue
#          as.integer(0), # exitstatus
#          as.double(0), # relerr
#          as.integer(0), # its
#	  PACKAGE="tweedie" )
#
#      density[i] <- tmp[[6]]
#		 its[i] <- tmp[[9]]
#
#}
#
#list( density=density, its=its)
#
#}



ptweedie_Series_Report <- function(q, power, mu, phi) {
  #
  # ptweedie.series
  #
  # ptweedie.series( y, power, mu, phi )
  #
  # Evaluating Tweedie cdf using a Poisson sum of incomplete gamma
  # functions.
  #
  # Arguments:
  #   y:  The values of a the random variable.  y  can be 
  #       a vector, and all elements must be non-negative.
  #
  #   mu: The (positive) mean of the distribution.
  #
  #   phi:    The (positive) dispersion parameter for the 
  #       distribution, such that the variance of the 
  #       observations is  var(Y) = phi * mu^power
  #
  #   power:  The index for the Tweedie distribution, such that 
  #       the variance of the observation is 
  #       var(Y) = phi * mu^power .  
  #
  # Description
  #   This function evaluates the Tweedie family of distributions
  #   distribution function by a Poisson sum of incomplete gammas functions.
  #
  # Value
  #   A list is returned with the following components:
  #   density:    The value of the density 
  #   
  # See Also:
  #   dtweedie.series.bigp (using series evaluation for evaluation)
  #   dtweedie (density evaluation using fast methods based on the
  #       series evaluation and cgf inversion)
  #
  # Examples
  #   ptweedie.series.smallp(mu=1, phi=1, y=data, power=3.3)
  
  # 
  # Peter K Dunn 
  # 14 March 2000 
  #
  y <- q
  # Error checks
  if ( power<1) stop("power must be between 1 and 2.\n")
  if ( power>2) stop("power must be between 1 and 2.\n")
  if ( any(phi<= 0)) stop("phi must be positive.\n")
  if ( any(y<0) ) stop("y must be a non-negative positive vector.\n")
  if ( any(mu<=0) ) stop("mu must be positive.\n")
  if ( length(mu)>1) {
    if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.\n")
  }else {
    mu <- array( dim=length(y), mu )
  }
  if ( length(phi)>1) {
    if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.\n")
  }else {
    phi <- array( dim=length(y), phi )
  }
  
  
  # Set up
  p <- power
  lambda <- mu^(2-p) / ( phi * (2-p) )
  tau <- phi * (p-1) * mu^( p-1 )
  alpha <- (2-p) / ( 1-p )
  
  
  # Now find the limits on N:
  # First the lower limit on N
  
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
  
  
  # Now for the upper limit on N
  
  lambda <- mu^(2-p) / ( phi * (2-p) )
  lambda <- min( lambda )
  logfmax <-  -log(lambda)/2
  estlogf <- logfmax
  
  N <- min( lambda )
  
  while ( estlogf > (logfmax - drop) ) {
    N <- N + 1
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  
  hi.N <- max( ceiling(N) )
  
  # Now evaluate between limits of N
  
  cdf <- array( dim=length(y), 0 )
  
  lambda <- mu^(2-p) / ( phi * (2-p) )
  tau <- phi * (p-1) * mu^( p-1 )
  alpha <- (2-p) / ( 1-p )
  
  
  for (N in (lo.N:hi.N)) {
    
    # Poisson density
    pois.den <- dpois( N, lambda)
    
    # Incomplete gamma
    incgamma.den <- pchisq(2*y/tau, -2*alpha*N )
    
    # What we want
    cdf <- cdf + pois.den * incgamma.den
    
  }
  
  cdf <- cdf + exp( -lambda )
  its <- hi.N - lo.N + 1
  
  
  list(cdf = cdf, 
       lo = lo.N, 
       hi = hi.N)
  
}





dtweedie_Inversion_Report <- function(y, power, mu, phi, exact=TRUE, method=3, verbose = FALSE){ 
  #
  # DTWEEDIE.INVERSION
  #
  # dtweedie.inversion( y, power, mu, phi, exact, rotate )
  #
  # Evaluating Tweedie densities using cgf inversion
  #
  # Arguments:
  #   y:  The values of a the random variable.  y  can be 
  #       a vector, and all elements must be non-negative.
  #
  #   mu: The (positive) mean of the distribution.
  #
  #   phi:    The (positive) dispersion parameter for the 
  #       distribution, such that the variance of the 
  #       observations is  var(Y) = phi * mu^power
  #
  #   power:  The index for the Tweedie distribution, such that 
  #       the variance of the observation is 
  #       var(Y) = phi * mu^power .  
  #
  #   exact:  If exact is TRUE, the exact zero acceleration algorithm 
  #       is used, otherwise the approx zero algorithm is used
  #
  #   method:  Which method to use.  If NULL, the optimal method is used
  #
  # Description
  #   This function evaluates the Tweedie family of distributions
  #   by inverting the cgf.
  #
  # Value
  #   A list is returned with the following components:
  #   density:    The value of the density 
  #   
  # See Also:
  #   dtweedie.series.bigp (using series evaluation for evaluation)
  #   dtweedie (density evaluation using fast methods based on the
  #       series evaluation and cgf inversion)
  #
  # Examples
  #   dtweedie.inversion(mu=1, phi=1, y=data, power=3.3)
  
  
  # 
  # Peter K Dunn 
  # 06 Aug 2002
  #
  
  # Error checks
  if ( power<1) stop("power must be greater than 1.")
  if ( any(phi<= 0)) stop("phi must be positive.")
  #if ( any(y<0) ) stop("y must be a non-negative vector.")
  if ( any(mu<=0) ) stop("mu must be positive.")
  if ( length(mu)>1) {
    if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
  }
  else {
    mu <- array( dim=length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi)>1) {
    if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
  }
  else {
    phi <- array( dim=length(y), phi )
    # A vector of all phi's
  }
  
  y.len <- length(y)
  density <- y
  its <- y
  
  if ( is.null(method)){
    method <- array( dim=length(y))
  } else {
    method <- array( method, dim=length(y))
  }
  
  for (i in (1:y.len)) {
    # Now choose the best method
    # There are three approaches, each a product of a simple bit
    # and a complicated bit computed in FORTRAN
    # We choose the method with the smaller 'simple' bit
    if ( y[i] <= 0 ) {
      if ( (power>1) & (power<2) ) {
        if ( y[i]==0 ) {
          density[i] <- exp( -mu[i] ^ (2-power) / ( phi[i] * (2-power) ) )
        } else {
          density[i] <- 0
        }
      } else {
        density[i] <- 0
      }
    } else {
      # Here, y > 0
      m2 <- 1/mu[i]
      
      theta <- ( mu[i]^(1-power) - 1 ) / ( 1 - power )
      if ( ( abs(power - 2 ) ) < 1.0e-07 ){
        kappa <- log(mu[i]) + (2 - power) * ( log(mu[i])^2 ) / 2
      } else {
        kappa <- ( mu[i]^(2-power) - 1 ) / ( 2 - power )
      }
      m1 <- exp( (y[i]*theta - kappa )/phi[i] )
      
      dev <- tweedie.dev(y=y[i], mu=mu[i], power=power )
      m3 <- exp( -dev/(2*phi[i]) ) / y[i]
      
      min.method <- min( m1, m2, m3 )
      
      # Now if no method requested, find the optimal      
      if ( is.null(method[i]) ) {
        if ( min.method==m1 ){
          use.method <- 1
        } else {
          if ( min.method==m2 ) {
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
      if ( verbose ) cat("Using method =",method,"\n")
      if ( use.method==2 ) {
        if ( verbose ) {
          cat("Sending y =",y[i]/mu[i],"; mu =",1,"; phi =",phi[i] / (mu[i]^(2-power)),"with method=",method[i],"\n")
        }
        tmp <- .Fortran( "twpdf",
                         as.double(power),
                         as.double(phi[i] / (mu[i]^(2-power)) ), # phi
                         as.double(y[i]/mu[i]), # y
                         as.double(1), # mu
                         as.integer( exact ), #exact as an integer
                         as.integer( verbose ), #verbose as an integer
                         as.double(0), # funvalue
                         as.integer(0), # exitstatus
                         as.double(0), # relerr
                         as.integer(0), # its
                         PACKAGE="tweedie")
        
        den <- tmp[[7]]
        density[i] <- den * m2
        method[i] <- 1
        its[i] <- tmp[[10]]
      } else {
        if ( use.method==1 ) {
          if ( verbose ) {
            cat("Sending y =",y[i],"; mu =",1,"; phi =",phi[i],"with method=",method[i],"\n")
          }
          tmp <- .Fortran( "pdf",
                           as.double(power),
                           as.double(phi[i]), # phi
                           as.double(y[i]), # y
                           as.double(1), # mu
                           as.integer( exact ), #exact as an integer
                           as.integer( verbose ), #verbose as an integer
                           as.double(0), # funvalue
                           as.integer(0), # exitstatus
                           as.double(0), # relerr
                           as.integer(0), # its
                           PACKAGE="tweedie")
          
          den <- tmp[[7]]
          density[i] <- den * m1
          method[i] <- 2
          its[i] <- tmp[[10]]
        } else { # use.method==3
          if ( verbose ) {
            cat("Sending y =",1,"; mu =",1,"; phi =",phi[i] / (y[i]^(2-power)),"with method=",method[i],"\n")
          }
          tmp <- .Fortran( "pdf",
                           as.double(power),
                           as.double(phi[i]/(y[i]^(2-power))), # phi
                           as.double(1), # y
                           as.double(1), # mu
                           as.integer( exact ), #exact as an integer
                           as.integer( verbose ), #verbose as an integer
                           as.double(0), # funvalue
                           as.integer(0), # exitstatus
                           as.double(0), # relerr
                           as.integer(0), # its
                           PACKAGE="tweedie")
          
          den <- tmp[[7]]
          density[i] <- den * m3
          method[i] <- 3
          its[i] <- tmp[[10]]
        }
      }
    }
  }
  
  list(density = density, 
       its = its, 
       method = method)
  
}

dtweedie_Inversion_Threemethods <- function(y, power, mu, phi, exact=TRUE){ 
  #
  # DTWEEDIE.INVERSION.THREEMETHODS
  #
  # dtweedie.inversion( y, power, mu, phi, exact, rotate )
  #
  # Evaluating Tweedie densities using cgf inversion
  #
  # Arguments:
  #   y:  The values of a the random variable.  y  can be 
  #       a vector, and all elements must be non-negative.
  #
  #   mu: The (positive) mean of the distribution.
  #
  #   phi:    The (positive) dispersion parameter for the 
  #       distribution, such that the variance of the 
  #       observations is  var(Y) = phi * mu^power
  #
  #   power:  The index for the Tweedie distribution, such that 
  #       the variance of the observation is 
  #       var(Y) = phi * mu^power .  
  #
  #   exact:  If exact is TRUE, the exact zero acceleration algorithm 
  #       is used, otherwise the approx zero algorithm is used
  #
  #   method:  Which method to use.  If NULL, the optimal method is used
  #
  # Description
  #   This function evaluates the Tweedie family of distributions
  #   by inverting the cgf.
  #
  # Value
  #   A list is returned with the following components:
  #   density:    The value of the density 
  #   
  # See Also:
  #   dtweedie.series.bigp (using series evaluation for evaluation)
  #   dtweedie (density evaluation using fast methods based on the
  #       series evaluation and cgf inversion)
  #
  # Examples
  #   dtweedie.inversion(mu=1, phi=1, y=data, power=3.3)
  
  
  # 
  # Peter K Dunn 
  # 06 Aug 2002
  #
  
  # Error checks
  if ( power<1) stop("power must be greater than 1.")
  if ( any(phi<= 0)) stop("phi must be positive.")
  if ( any(y<0) ) stop("y must be a non-negative vector.")
  if ( any(mu<=0) ) stop("mu must be positive.")
  if ( length(mu)>1) {
    if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
  }
  else {
    mu <- array( dim=length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi)>1) {
    if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
  }
  else {
    phi <- array( dim=length(y), phi )
    # A vector of all phi's
  }
  
  y.len <- length(y)
  density1 <- density2 <- density3 <- y
  its1 <- its2 <- its3 <- y
  
  
  for (i in (1:y.len)) {
    
    # Now choose the best method
    # There are three approaches, each a product of a simple bit
    # and a complicated bit computed in FORTRAN
    # We choose the method with the smaller 'simple' bit
    
    m1 <- 1/mu[i]
    
    theta <- ( mu[i]^(1-power) - 1 ) / ( 1 - power )
    kappa <- ( mu[i]^(2-power) - 1 ) / ( 2 - power )
    m2 <- exp( (y[i]*theta - kappa )/phi[i] )
    
    dev <- tweedie.dev(y=y[i], mu=mu[i], power=power )
    m3 <- exp( -dev/(2*phi[i]) ) / y[i]
    
    # Now use the methods
    # NOTE: FOR ALL  METHODS, WE HAVE mu=1
    
    tmp1 <- .Fortran( "pdf",
                      as.double(power),
                      as.double(phi[i] / (mu[i]^(2-power)) ), # phi
                      as.double(y[i]/mu[i]), # y
                      as.double(1), # mu
                      as.integer( exact ), #exact as an integer
                      as.double(0), # funvalue
                      as.integer(0), # exitstatus
                      as.double(0), # relerr
                      as.integer(0), # its
                      PACKAGE="tweedie")
    
    den <- tmp1[[6]]
    density1[i] <- den * m1
    method <- 1
    its1[i] <- tmp1[[9]]
    
    tmp2 <- .Fortran( "pdf",
                      as.double(power),
                      as.double(phi[i]), # phi
                      as.double(y[i]), # y
                      as.double(1), # mu
                      as.integer( exact ), #exact as an integer
                      as.double(0), # funvalue
                      as.integer(0), # exitstatus
                      as.double(0), # relerr
                      as.integer(0), # its
                      PACKAGE="tweedie")
    
    den <- tmp2[[6]]
    density2[i] <- den * m2
    its2[i] <- tmp2[[9]]
    
    tmp3 <- .Fortran( "pdf",
                      as.double(power),
                      as.double(phi[i]/(y[i]^(2-power))), # phi
                      as.double(1), # y
                      as.double(1), # mu
                      as.integer( exact ), #exact as an integer
                      as.double(0), # funvalue
                      as.integer(0), # exitstatus
                      as.double(0), # relerr
                      as.integer(0), # its
                      PACKAGE="tweedie")
    
    den <- tmp3[[6]]
    density3[i] <- den * m3
    its3[i] <- tmp3[[9]]
  }
  list(density1 = density1, 
       density2 = density2, 
       density3 = density3,
       m1 = m1, 
       m2 = m2, 
       m3 = m3,
       fortran1 = tmp1[[6]], 
       fortran2 = tmp2[[6]], 
       fortran3 = tmp3[[6]],
       its1 = its1, 
       its2 = its2, 
       its3 = its3)
}


dtweedie_Fortran <- function(y, power, phi, exact=TRUE){ 
  #
  # DTWEEDIE.FORTRAN
  #
  # dtweedie.inversion( y, power, mu, phi, exact, rotate )
  #
  # Evaluates the FORTRAN directly, and always uses mu=1
  #
  # Arguments:
  #   y:  The values of a the random variable.  y  can be 
  #       a vector, and all elements must be non-negative.
  #
  #   mu: The (positive) mean of the distribution.
  #
  #   phi:    The (positive) dispersion parameter for the 
  #       distribution, such that the variance of the 
  #       observations is  var(Y) = phi * mu^power
  #
  #   power:  The index for the Tweedie distribution, such that 
  #       the variance of the observation is 
  #       var(Y) = phi * mu^power .  
  #
  #   exact:  If exact is TRUE, the exact zero acceleration algorithm 
  #       is used, otherwise the approx zero algorithm is used
  #
  #   method:  Which method to use.  If NULL, the optimal method is used
  #
  # Description
  #   This function evaluates the Tweedie family of distributions
  #   by inverting the cgf.
  #
  # Value
  #   A list is returned with the following components:
  #   density:    The value of the density 
  #   
  # See Also:
  #   dtweedie.series.bigp (using series evaluation for evaluation)
  #   dtweedie (density evaluation using fast methods based on the
  #       series evaluation and cgf inversion)
  #
  # Examples
  #   dtweedie.inversion(mu=1, phi=1, y=data, power=3.3)
  
  
  # 
  # Peter K Dunn 
  # 06 Aug 2002
  #
  
  # Error checks
  if ( power<1) stop("power must be greater than 1.")
  if ( any(phi<= 0)) stop("phi must be positive.")
  if ( any(y<0) ) stop("y must be a non-negative vector.")
  mu <- array( dim=length(y), 1 )
  # A vector of all mu's (all one)
  
  if ( length(phi)>1) {
    if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
  }
  else {
    phi <- array( dim=length(y), phi )
    # A vector of all phi's
  }
  
  y.len <- length(y)
  density <- y
  its <- y
  
  
  for (i in (1:y.len)) {
    
    # NOTE: FOR ALL  METHODS, WE HAVE mu=1
    tmp <- .Fortran( "pdf",
                     as.double(power),
                     as.double(phi[i] ), # phi
                     as.double(y[i]), # y
                     as.double(1), # mu
                     as.integer( exact ), #exact as an integer
                     as.double(0), # funvalue
                     as.integer(0), # exitstatus
                     as.double(0), # relerr
                     as.integer(0), # its
                     PACKAGE="tweedie")
    
    den <- tmp[[6]]
    density[i] <- den
    its[i] <- tmp[[9]]
  }
  list(density = density, 
       its = its)
  
}




# 
# 
# 
# 
# dtweedie.dldphi <- function(phi, mu, power, y ){
# #
# # DTWEEDIE.DLDPHI
# #
# # This S-Plus function calculates the log-likelihood
# # function, wrt phi, for p>2.  In particular, it returns
# #    sum{ d(log f)/d(phi) } = d( log-likelihood )/d(phi).
# # The mle of phi can be found, therefore, by setting this to zero.
# # y  is generally a vector of observed values.
# 
# #
# # Peter Dunn
# # 31 Jan 2001
# #
# 
# if ( (power != 2 ) & ( power != 1 ) ) {
# 
#    k <- phi^(1/(power-2))
# 
#    if ( k < 1 ) {
#       # Use the transform f(y; mu, phi) = c f(c*y; c*mu, c^(2-p)*phi)
#       # and differentiate with c=phi^(1/(p-2)):
#       #    d log f / d phi = c^(2-p) * {df(cy; c*mu, 1)/dphi} / f(cy; c*mu, 1)
#       f <- dtweedie( y=k*y, power=power, mu=k*mu, phi=1 )
#       d <- dtweedie.dlogfdphi( y=k*y, power=power, mu=k*mu, phi=1 )
#          # Note:  We need dlogf/dphi = dlogf.dphi * f
#       top <- d * f
#       d <- -2* sum( top / f * k^(2-power) )
# 
#    }
#    else{
#     # Compute directly
#       d <- -2*sum( dtweedie.dlogfdphi(y=y, power=power, mu=mu, phi=phi) )
#    }
# }
# else{
#    # Cases p==1 and  p=2
#    d <- -2*sum( dtweedie.dlogfdphi(y=y, power=power, mu=mu, phi=phi) )
# }
# d
# }

tweedie_Igrand <- function(phi, power, y, mu = ifelse(type=="PDF", 1, NULL), t, doPlot = TRUE, type = "PDF" ){
  
  if ( type == "PDF"){
    # DENSITY FUNCTION
    # Assumes mu=1

    front <- 1/( phi * ( 2 - power ) )
    alpha <- (2 - power)/(1 - power)
    omega <- atan( (1 - power) * t * phi)

    k.real <- front * cos(omega * alpha)/( cos(omega)^alpha ) - front
    k.imag <- front * sin(omega * alpha)/( cos(omega)^alpha ) - y * t
    
    if ( power > 2 ) {
      ig <- exp( k.real) * cos( k.imag )
    }
    if ( (power>1) & (power < 2 ) ) {
      lambda1 <- (1 / (phi * (2 - power)))
      ig <- exp( k.real) * cos( k.imag ) - exp(-lambda1) * cos(t * y)
    }
  } else {
    # DISTRIBUTION FUNCTION   
    front <- mu ^ (2 - power) /( phi * ( 2 - power ) )
    alpha <- (2 - power)/(1 - power)
    omega <- atan( (1 - power) * t * phi / mu^(1 - power))
    
    k.real <- front * (cos(omega * alpha)/( cos(omega)^alpha ) - 1) 
    k.imag <- front * sin(omega * alpha)/( cos(omega)^alpha ) - y * t
    
    if ( power > 2 ) {
      ig <- exp( k.real) * cos( k.imag )
    }
    if ( (power > 1) & (power < 2 ) ) {
      lambda1 <- (1 / (phi * (2 - power)))
      ig <- exp( k.real) * cos( k.imag ) / t
    }
  }
  if ( doPlot){
    if (type == "PDF") {
      plot( y = exp(k.real) * cos(k.imag),
            x = t,
            lwd = 2,
            las = 1,
            type = "l", 
            main = expression( atop( bold(Integrand~PDF),
                                     ("for"~mu==1))),
            xlab = expression(italic(t)),
            ylab = "Integrand")
      abline(h = 0,
             col = "grey")
      lines( y = exp(k.real),
             x = t,
             lwd = 1,
             col = "grey")
      lines( y = -exp(k.real),
             x = t,
             lwd = 1,
             col = "grey")      
    } else {
      plot( y = exp(k.real) * sin(k.imag) / t,
            x = t,
            lwd = 2,
            las = 1,
            type = "l", 
            main = paste("Integrand for DF"),
            xlab = expression(italic(t)),
            ylab = "Integrand")
      abline(h = 0,
             col = "grey")
      lines( y = exp(k.real)/t,
             x = t,
             lwd = 1,
             col = "grey")
      lines( y = -exp(k.real)/t,
             x = t,
             lwd = 1,
             col = "grey")
    }    
  }
  
  list(ig = ig, 
       k.real = k.real, 
       k.imag = k.imag,
       type = type)
}
