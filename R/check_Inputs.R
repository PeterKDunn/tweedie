check_Inputs <- function(y, mu, phi, power){
  
  y.negative <- (y < 0)
  y.original <- y
  y <- y[ !y.negative ]

  inputs_Error <- FALSE
  
  # Checking the input values
  if ( any(power < 1) )  {
    stop( "The Tweedie index parameter must be greater than 1.\n")
    inputs_Error <- TRUE
  }
  if ( any(phi <= 0) ) {
    stop("phi must be positive.")
    inputs_Error <- TRUE
  }
  if ( any(mu <= 0) ) {
    stop("mu must be positive.\n")
    inputs_Error <- TRUE
  }
  
  # Checking the input length of  mu
  if ( length(mu) > 1) {
    # If  mu  not a scalar, check it is the same length as  y

    if ( length(mu) != length(y.original) ) stop("mu must be scalar, or the same length as y.\n")
    inputs_Error <- TRUE
    
  } else {
    # If  mu  is a scalar, force it to be the same same length as  y
    
    mu <- array( dim = length(y.original), mu )
    # A vector of all mu's
  
  }
  

  # Checking the input length of  phi  and  mu
  # The length of xi/power should have been checked in sort_Notation, so does not need checking here again
  
  if ( length(phi) > 1) {
    # If  phi  not a scalar, check it is the same length as  y

    if ( length(phi) != length(y.original) ) stop("phi must be scalar, or the same length as y.\n")
    inputs_Error <- TRUE
    
  } else {
    # If  mu  is a scalar, force it to be the same same length as  y

    phi <- array( dim = length(y.original), phi )
    # A vector of all phi's
    
  }
  
  
  # Now remove values where  y  is negative
  mu.original <- mu
  mu <- mu[ !y.negative ]
  
  phi.original <- phi
  phi <- phi[ !y.negative ]
  
  return( list(mu.original = mu.original,
               mu = mu, 
               phi.original = phi.original,
               phi = phi,
               inputs_Error = inputs_Error) )
} 