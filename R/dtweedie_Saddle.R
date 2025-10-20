dtweedie.saddle <- function(y, xi=NULL, mu, phi, eps=1/6, power=NULL) {
  #
  # Peter K Dunn
  # 09 Jan 2002
  #
  #
  
  # Sort out the xi/power notation
  if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
  xi.notation <- TRUE
  if ( is.null(power) ) {
    power <- xi
  } else {
    xi.notation <- FALSE
  }
  if ( is.null(xi) ) {
    xi.notation <- FALSE
    xi <- power
  }
  if ( xi != power ) {
    cat("Different values for xi and power given; the value of xi used.\n")
    power <- xi
  }
  index.par       <- ifelse( xi.notation, "xi", "p")
  index.par.long  <- ifelse( xi.notation, "xi", "power")
  
  
  # Error traps
  #
  if( any(phi <= 0) )
    stop("phi must be positive.")
  if( (power >= 1) & any(y < 0))
    stop("y must be a non-negative vector.")
  if( any(mu <= 0) )
    stop("mu must be positive.")
  if ( length(mu) == 1 )
    mu <- array( dim=length(y), mu )
  if ( length(phi) == 1 )
    phi <- array( dim=length(y), phi )
  #
  # END error checks
  #
  #   if ( any(y==0) ) y <- y + (1/6)
  y.eps <- y 
  
  if (power<2) y.eps <- y + eps
  # Nelder and Pregibon's idea for problems at y=0
  
  y0 <- (y == 0)
  density <- y
  #    if(any(y == 0)) {
  #        if(power >= 2) {
  #         # This is a problem!
  #            density[y0] <- 0 * y[y0]
  #        }
  #      else{
  
  dev <- tweedie.dev(y = y, 
                     mu = mu, 
                     power = power)
  density <- (2 * pi * phi * y.eps ^ power) ^ (-1/2) * exp( -dev / (2 * phi) )
  
  #      }
  #   }
  if ( any(y == 0) ){
    if((power >= 1) && (power < 2)) {
      lambda <- mu[y0] ^ (2 - power) / (phi[y0] * (2 - power))
      density[y0] <- exp( -lambda)
    } else {
      density[y0] <- 0
    }
  }
  
  #    if(any(y != 0)) {
  #     dev <- tweedie.dev(y=y[yp], mu=mu[yp], power=power)
  #      density[yp] <- (2*pi*phi[yp]*y[yp]^power)^(-1/2) * exp( -dev/(2*phi[yp]) )
  #    }
  density
}


#############################################################################
