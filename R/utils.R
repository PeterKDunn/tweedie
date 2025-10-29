sort_Notation <- function(xi = NULL, power = NULL){
  # Sorts out whether the  xi  or the  p/power  notation is being used,
  # so output presented appropriately
  
  if ( length(xi) > 1) stop(" The Tweedie index parameter (xi) must be a single value.")
  if ( length(power) > 1) stop(" The Tweedie index paramerer (power) must be a single value.")
  
  if ( is.null(xi) & is.null(power) ) stop("The Tweedie index parameter (xi) must be given.")
  
  if ( !is.null(xi) & !is.null(power) ){
    # That is: if BOTH xi and power are given
    if (xi == power) {
      power <- NULL # Either is OK, but default to  xi  notation ...and carry it through below
    } else {
      if (xi != power) {
        power <- NULL # Use the value if xi, and carry of through below
        message("Both  xi  and  power  are given; using the given value of  xi")
      }
    }
  }
  
  if ( is.null(xi) & !is.null(power)) {
    # xi  is NULL
    xi <- power
    xi.notation <- FALSE
    index.par <- "p"
    index.par.long <- "power"
  } 
  if ( !is.null(xi) & is.null(power)) {
    # power  is NULL
    power <- xi
    xi.notation <- TRUE
    index.par <- "xi"
    index.par.long <- "xi"
  }
  
  return( list(xi = xi, 
               power = power, 
               xi.notation = xi.notation, 
               index.par = index.par, 
               index.par.long = index.par.long) )
} 



################################################################################

check_Inputs <- function(y, mu, phi, power){
  # Checks that the inputs satisfy the necessary criteria (e.g., mu > 0).
  # Ensures that y, mu and phi are all vectors of the same length.
  
  y.original <- y
  y.negative <- (y < 0)
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