dtweedie <- function(y, xi = NULL, mu, phi, power = NULL, verbose = FALSE, details = FALSE){

  # Evaluates the density for Tweedie distributions, for given values of:
  #   y (possibly a vector)
  #   mu, the mean 
  #   phi the dispersion parameter
  #   power,  the Tweedie index parameter
  #   verbose: whether to display what is happening
  #   details:  whether to returns reports of relerr, regions needed, etc.

  # Two methods are employed:  cgf inversion (type=1)
  # and series evaluation (type=2 if 1<p<2; type=3 if p>2)).
  #
  # This function uses bivariate interpolation to accurately
  # approximate the inversion in conjunction with the series.
  #
  #
  # FIRST, establish whether the series or the cgf inversion is the better method.
  
  #
  # Here is a summary of what happens:
  #
  #   p                          |  What happpens
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
  
  ### BEGIN preliminary work
  
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
  out <- check_Inputs(y, mu, phi, power)
  mu <- out$mu
  phi <- out$phi
  f <- array(0,
             dim = length(q) )
  if (details) regions <- array(0, dim = length(q))
  
  
  # IDENTIFY SPECIAL CASES
  special_y_Cases <- rep(FALSE, length(y))
  if (verbose) cat("- Checking for special cases\n")
  out <- special_Cases(y, mu, phi, power)
  special_p_Cases <- out$special_p_Cases
  special_y_Cases <- out$special_y_Cases
  if (verbose & special_p_Cases) cat("  - Special case for p used\n")
  if ( any(special_y_Cases) ) {
    special_y_Cases <- out$special_y_Cases  
    if (verbose) cat("  - Special cases for first input found\n")
    f <- out$f
  }
  
  ### END preliminary work
  
  
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
