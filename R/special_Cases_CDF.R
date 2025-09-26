special_Cases_CDF <- function(y, mu, phi, power){
  
  f <- array( dim = length(y) )

  # Special Cases
  if ( power == 2 ) {
    f <- pgamma( rate = 1 / (phi * mu), 
                 shape = 1 / phi, 
                 q = y )
  }
  if ( power == 0) {
    f <- pnorm( mean = mu, 
                sd = sqrt(phi),
                q = y )
  }
  if ( power == 1) {
    f <- ppois(q = y, 
               lambda = mu / phi )
  }
  return( list(f = f) )
}
