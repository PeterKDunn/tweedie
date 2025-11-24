

#############################################################################
dtweedie_stable <- function(y, power, mu, phi)
{
  # Error checks
  if ( power < 1) stop("power must be greater than 2.\n")
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
  alpha <- (2 - power) / (1 - power)
  beta <- 1
  k <- 1  # The parameterization used
  delta <- 0
  gamma <- phi * (power - 1) * 
    ( 1 / (phi * (power - 2)) * cos( alpha * pi / 2 ) ) ^ (1 / alpha)
  
  #        require(stabledist) # Needed for  dstable
  ds <- stabledist::dstable(y,
                            alpha = alpha, 
                            beta = beta, 
                            gamma = gamma, 
                            delta = delta, 
                            pm = k)
  density <- exp((y * mu ^ (1 - power) / (1 - power) - mu ^ (2 - power) / (2 - power)) / phi) * ds
  
  
  density
}




#############################################################################
tweedie_AIC <- function( glm.obj, dispersion = NULL, k = 2, verbose = TRUE){ 
  # New  dispersion  input for (e.g.) Poisson case, added October 2017
  
  wt <- glm.obj$prior.weights
  n <- length(glm.obj$residuals)
  edf <- glm.obj$rank  # As used in logLik.glm()
  
  mu <- stats::fitted( glm.obj )
  y  <- glm.obj$y
  p <- get("p", envir = environment(glm.obj$family$variance))
  
  if ( is.null(dispersion)) {  # New section
    if (p == 1 & verbose) message("*** Tweedie index power = 1: Consider using  dispersion=1  in call to  AICtweedie().\n")
    dev <- deviance(glm.obj)
    disp <- dev / sum(wt)  # In line with Gamma()$aic
    edf <- edf + 1  # ADD one as we are estimating phi too
  } else {
    disp <- dispersion
  }
  
  den <- dtweedie( y = y, 
                   mu = mu, 
                   phi = disp, 
                   power = p)
  AIC <- -2 * sum( log(den) * wt) 
  
  return( AIC + k * (edf) )
  
}




#' @export
AICtweedie <- function( glm.obj, dispersion = NULL, k = 2, verbose = TRUE){ 
  .Deprecated("tweedie_AIC", package = "tweedie")
  tweedie_AIC(glm.obj = glm.obj, 
              dispersion = dispersion,
              k = k,
              verbose = verbose)
}



#############################################################################

tweedie_convert <- function(xi = NULL, mu, phi, power = NULL){
  ### ADDED 14 July 2017
  
  if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
  xi.notation <- TRUE
  if ( is.null(power) ) {   # Then  xi  is given
    if ( !is.numeric(xi)) stop("xi  must be numeric.\n")
    power <- xi
  } else {
    xi.notation <- FALSE
  }
  if ( is.null(xi) ) {   # Then   power  is given
    if ( !is.numeric(power)) stop("power  must be numeric.\n")
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
  if ( power < 1)      stop( paste(index.par.long, "must be greater than 1.\n") )
  if ( power >= 2)     stop( paste(index.par.long, "must be less than 2.\n") )
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(mu <= 0) )  stop("mu must be positive.\n")
  
  if( length(mu) != length(phi) ){
    if ( length(mu)  == 1 ) mu  <- array(dim = length(phi), mu  ) 
    if ( length(phi) == 1 ) phi <- array(dim = length(mu),  phi ) 
  }
  # Now  mu  and  phi  will be the same length if one of them happened to be a scalar, so this works:
  if( length(mu) != length(phi) ) stop("phi and mu must be scalars, or the same length.\n")
  
  lambda <- ( mu ^ (2 - xi) ) / ( phi * (2 - xi) )  # Poisson distribution mean
  alpha  <- (2 - xi)  / (xi - 1)             # gamma distribution alpha (shape)
  gam    <- phi * (xi - 1) * mu ^ (xi - 1)   # gamma distribution beta  (scale)
  p0     <- exp( -lambda )
  phi.g  <- (2 - xi) * (xi - 1) * phi ^ 2 * mu ^ ( 2 * (xi-1) )
  mu     <- gam / phi
  
  list( poisson.lambda = lambda, 
        gamma.shape = alpha, 
        gamma.scale = gam, 
        p0 = p0,
        gamma.mean = mu,
        gamma.phi = phi.g)
  
}




#' @export
tweedie.convert <- function(xi = NULL, mu, phi, power = NULL){
  .Deprecated("tweedie_convert", package = "tweedie")
  if (is.null(power)) power <- xi
  tweedie_convert(xi = NULL, 
              mu = mu,
              phi = phi,
              power=power)
}
