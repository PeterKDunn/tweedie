kappa <- function(mu, p){
  if (p == 2) {
    k <- (mu^(2 - p) - 1) / (1 - p)  
  } else {
    k <- log(mu)    
  }
  return(k)
}



theta <- function(mu, p){
  if (p == 1) {
    theta <- (mu^(1 - p) - 1) / (1 - p)  
  } else {
    theta <- log(mu)    
  }
  return(theta)
}




K <- function(theta, phi, t){
  (kappa(theta + t*phi) - kappa(theta)) / phi
}



k <- function(p, mu, phi, y, t, verbose=FALSE){
  front <- ( mu^(2 - p) ) / (phi * (2 - p) )
  alpha <- (2 - p)/(1 - p) 
  omega <- atan( ( (1 - p) * t * phi) / mu^(1 - p) )

  if (verbose){
    cat("omega: ", omega, "\n")
  }
  if ( any( (omega <= 0) & (omega >= -pi/2) ) ) {
    # ALL OK
  } else {
    cat("OOPS: change atan() call??\n")
  }

  Real <- front * ( cos(omega * alpha) / ( cos(omega)^alpha) - 1) 
  Imag <- front *   sin(omega * alpha) / ( cos(omega)^alpha ) - (t * y)
  
  return( list(Real = Real,
               Imag = Imag) )
}


kdash <- function(p, mu, phi, y, t){
  omega <- atan( ( (1 - p) * t * phi) / mu^(1 - p) )
  
  if ( any( (omega <= 0) & (omega >= -pi/2) ) ) {
    # ALL OK
  } else {
    cat("OOPS: change atan() call??\n")
  }
  
  mu * cos(omega / (1 - p) ) / (cos(omega)^(1/(1 - p))) - y

}


kdashdash <- function(p, mu, phi, y, t){
  omega <- atan( ( (1 - p) * t * phi) / mu^(1 - p) )
  
  -phi * mu ^ (p) * sin(omega * p / (1 - p) ) / (cos(omega)^(p/(1 - p))) 
  
}

kdashdashdash <- function(p, mu, phi, y, t){
  omega <- atan( ( (1 - p) * t * phi) / mu^(1 - p) )
  
  -phi^2 * p  * mu^(2*p - 1) * cos(omega * (2*p - 1) / (1 - p) ) / (cos(omega)^((2*p - 1)/(1 - p))) 
  
}



igrand <- function(p, mu, phi, y, t){
  rk <- k(p = p,
          mu = mu, 
          phi = phi, 
          y = y, 
          t = t)
  cat("t, Rek, Imk", t, rk$Real, rk$Imag, "\n")
  igrand <- exp( rk$Real ) * sin(rk$Imag) / t
  if ( p < 2 ) {
    igrand <- igrand - ( exp(rk$Real) * sin(rk$Imag - t*y) ) / t
  }
  return(igrand)
}





igrandPDF <- function(p, mu, phi, y, t){
  rk <- k(p = p,
          mu = mu, 
          phi = phi, 
          y = y, 
          t = t)
#  cat("t, Rek, Imk", t, rk$Real, rk$Imag, "\n")
  igrand <- exp( rk$Real ) * cos(rk$Imag)
#  if ( p < 2 ) {
#    igrand <- igrand - ( exp(rk$Real) * sin(rk$Imag - t*y) ) / t
#  }
  return(igrand)
}

