
#############################################################################
rtweedie <- function(n, xi = NULL, mu, phi, power = NULL){
  
  out <- sort_Notation(xi = NULL, power = NULL)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long ### MAY NOT BE NEEDED!!!
  

  # Error checks
  if ( any(power < 1) )  stop( paste(index.par.long, "must be greater than 1.\n") )
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( n < 1 ) stop("n must be a positive integer.\n")
  if ( any(mu <= 0) ) stop("mu must be positive.\n")
  if ( length(mu) > 1) {
    if ( length(mu) != n ) stop("mu must be scalar, or of length n.\n")
  } else {
    mu <- array( dim = n, mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != n ) stop("phi must be scalar, or of length n.\n")
  } else {
    phi <- array( dim = n, phi )
    # A vector of all phi's
  }
  if (power==1) {
    rt <- phi * rpois(n, 
                      lambda = mu / ( phi  ) )
  }
  
  if (power == 2) {
    alpha <- (2 - power) / (1 - power)
    gam <- phi * (power - 1) * mu ^ (power - 1)
    rt <- rgamma( n, 
                  shape = 1 / phi, 
                  scale = gam )
  }
  
  if ( power > 2) {
    rt <- qtweedie( runif(n),
                    mu = mu,
                    phi = phi, 
                    power = power)
  }
  
  if ( (power > 1) & (power < 2) ) {
    # Two options:  As above or directly.
    # Directly is faster
    rt <- array( dim = n, NA)
    
    lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
    alpha <- (2 - power) / (1 - power)
    gam <- phi * (power - 1) * mu ^ (power - 1)
    
    N <- rpois(n, 
               lambda = lambda)
    for (i in (1:n) ){
      rt[i] <- rgamma(1, 
                      shape = -N[i] * alpha, 
                      scale = gam[i])
    }
  }
  as.vector(rt)
}
