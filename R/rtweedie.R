rtweedie <- function(n, xi = NULL, mu, phi, power = NULL){
  
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
  out <- check_Inputs(q, mu, phi, power,
                      random_Numbers = TRUE)
  mu <- out$mu
  phi <- out$phi
  f <- array(0,
             dim = length(q) )
  if (details) regions <- array(0, dim = length(q))
  
  
  # IDENTIFY SPECIAL CASES
  if (verbose) cat("- Checking for special cases\n")
  out <- special_Cases(q, mu, phi, power)
  f <- out$f
  special_p_Cases <- out$special_p_Cases
  if (verbose & special_p_Cases) cat("  - Special case for p used\n")
  special_y_Cases <- out$special_y_Cases
  if (verbose & any(special_y_Cases)) cat("  - Special cases for first input found\n")
  
  ### END preliminary work
  
  
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
