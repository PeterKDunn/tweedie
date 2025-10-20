dtweedie <- function(y, xi = NULL, mu, phi, power = NULL, verbose = FALSE, details = FALSE)
{
  #
  # This is a function for determining Tweedie densities.
  # Two methods are employed:  cgf inversion (type=1)
  # and series evaluation (type=2 if 1<p<2; type=3 if p>2)).
  # This function uses bivariate interpolation to accurately
  # approximate the inversion in conjunction with the series.
  #
  #
  # FIRST, establish whether the series or the cgf
  # # inversion is the better method.
  #
  
  #
  # Here is a summary of what happens:
  #
  #   p                          |  Whats happpens
  # -----------------------------+-------------------
  #  1 -> p.lowest (=1.3)        |  Use series A
  #  p.lowest (=1.3) -> p=2      |  Use series A if
  #                              |        xix > xix.smallp (=0.8)
  #                              |      inversion with rho if
  #                              |        xix < xix.smallp (=0.8)
  #  p=2 -> p=3                  |  Use series B if
  #                              |        xix > xix.smallp (=0.8)
  #                              |      inversion with rho if
  #                              |        xix < xix.smallp (=0.8)
  #  p=3 -> p.highest (=4)       |  Use series B if xix > xix.bigp (=0.9)
  #                              |      inversion with 1/rho if
  #                              |        xix < xix.bigp (=0.9)
  #  p > 4                       |  Use series B if xix > 0.03
  #                              | Saddlepoint approximation if xix<= 0.03 (NEEDS FIXING)
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
  
  
  # Error checks
  if ( any(power < 1) )  stop( paste(index.par.long, "must be greater than 1.\n") )
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y < 0) ) stop("y must be a non-negative vector.\n")
  if ( any(mu <= 0) ) stop("mu must be positive.\n")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.\n")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.\n")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  density <- y
  
  # Special Cases
  if ( power == 3 ){
    density <- statmod::dinvgauss(x = y, 
                                  mean = mu, 
                                  dispersion = phi)
    return(density)
  }
  if ( power == 2 ) {
    density <- dgamma( rate = 1 / (phi * mu), 
                       shape = 1 / phi, 
                       x = y )
    return(density)
  }
  if ( power == 0) {
    density <- dnorm( mean = mu, 
                      sd = sqrt(phi), 
                      x = y )
    return(density)
  }
  if ( (power == 1) & (all(phi == 1))) {
    # Poisson case
    density <- dpois(x = y / phi, 
                     lambda = mu / phi )
    return(density)
  }
  
  # Set up
  id.type0 <- array( dim = length(y) )
  id.series <- id.type0
  id.interp <- id.type0
  
  ###   Now consider the cases 1 < p < 2   ###
  id.type0 <- (y == 0)
  if (any(id.type0) ) {
    if (power > 2) {
      density[id.type0] <- 0
    } else {
      lambda <- mu[id.type0] ^ (2 - power) / (phi[id.type0] * (2 - power))
      density[id.type0] <- exp( -lambda )
    }
  }
  
  
  xi <- array( dim = length(y) )
  xi[id.type0] <- 0
  xi[!id.type0] <- phi[!id.type0] * y[!id.type0] ^ (power - 2)
  xix <- xi / ( 1 + xi )
  
  if ( (power > 1) && (power <= 1.1) ) {
    id.series <- (!id.type0)
    if (any(id.series)){
      density[id.series] <- dtweedie.series(y = y[id.series],
                                            mu = mu[id.series], 
                                            phi = phi[id.series],
                                            power = power)
    }
    return(density = density)
  }
  
  if ( power == 1 ) { # AND phi not equal to one here
    id.series <- rep(TRUE, length(id.series))
    id.interp <- rep(FALSE, length(id.series))
  }
  if ( (power > 1.1) && (power <= 1.2) ) {
    id.interp <- ( (xix > 0) & (xix < 0.1) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.1
      p.hi <- 1.2
      xix.lo <- 0
      xix.hi <- 0.1
      np <- 15
      nx <- 25
    }
  }
  
  if ( (power > 1.2) && (power <= 1.3) ) {
    id.interp <- ( (xix > 0) & (xix < 0.3) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.2
      p.hi <- 1.3
      xix.lo <- 0
      xix.hi <- 0.3
      np <- 15
      nx <- 25
    }
  }
  if ( (power > 1.3) && (power <= 1.4) ) {
    id.interp <- ( (xix > 0) & (xix < 0.5) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.3
      p.hi <- 1.4
      xix.lo <- 0
      xix.hi <- 0.5
      np <- 15
      nx <- 25
    }
  }
  if ( (power > 1.4) && (power <= 1.5) ) {
    id.interp <- ( (xix > 0) & (xix < 0.8) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.4
      p.hi <- 1.5
      xix.lo <- 0
      xix.hi <- 0.8
      np <- 15
      nx <- 25
    }
  }
  if ( (power > 1.5) && (power < 2) ) {
    id.interp <- ( (xix > 0) & (xix < 0.9) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.5
      p.hi <- 2
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
    }
  }
  
  ### Cases p>2   ###
  if ( (power > 2) && (power < 3) ) {
    id.interp <- ( (xix > 0) & (xix < 0.9) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 2
      p.hi <- 3
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
    }
  }
  if ( (power >= 3) && (power < 4) ) {
    id.interp <- ( (xix > 0) & (xix < 0.9) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 3
      p.hi <- 4
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
    }
  }
  if ( (power >= 4) && (power < 5) ) {
    id.interp <- ( (xix > 0) & (xix < 0.9) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 4
      p.hi <- 5
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
    }
  }
  if ( (power >= 5) && (power < 7) ) {
    id.interp <- ( (xix > 0) & (xix < 0.5) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 5
      p.hi <- 7
      xix.lo <- 0
      xix.hi <- 0.5
      np <- 15
      nx <- 25
    }
  }
  if ( (power >= 7) && (power <= 10) ) {
    id.interp <- ( (xix > 0) & (xix < 0.3) )
    id.series <- (!(id.interp | id.type0))
    if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 7
      p.hi <- 10
      xix.lo <- 0
      xix.hi <- 0.3
      np <- 15
      nx <- 25
    }
  }
  
  if ( power > 10) {
    id.series <- (y != 0)
    id.interp <- (!(id.series | id.type0))
  }
  
  if (any(id.series)) {
    density[id.series] <- dtweedie.series(y = y[id.series],
                                          mu = mu[id.series], 
                                          phi = phi[id.series],
                                          power = power)
  }
  
  if (any(id.interp)) {
    dim( grid ) <- c( nx + 1, np + 1 )
    rho <- dtweedie.interp(grid, 
                           np = np, 
                           nx = nx,
                           xix.lo = xix.lo, 
                           xix.hi = xix.hi,
                           p.lo = p.lo, 
                           p.hi = p.hi,
                           power = power, 
                           xix = xix[id.interp] )
    dev <- tweedie.dev(power = power, 
                       mu = mu[id.interp], 
                       y = y[id.interp])
    front <- rho / (y[id.interp] * sqrt(2 * pi * xi[id.interp]))
    density[id.interp] <- front * exp(-1/(2 * phi[id.interp]) * dev)
  }
  
  
  density
  
}
