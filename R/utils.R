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
  
  inputs_Error <- FALSE
  
  # Checking the input values: power
  if ( any(power < 1) )  {
    stop( "The Tweedie index parameter must be greater than 1.\n")
    inputs_Error <- TRUE
  }
  
  # Checking the input values: phi
  if ( any(phi <= 0) ) {
    stop("phi must be positive.")
    inputs_Error <- TRUE
  }
  
  # Checking the input values: mu
  if ( any(mu <= 0) ) {
    stop("mu must be positive.\n")
    inputs_Error <- TRUE
  }
  
  
  # Checking the length of  mu
  if ( length(mu) > 1) {
    # If  mu  not a scalar, check it is the same length as  y
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as the first input.\n")
    inputs_Error <- TRUE
    
  } else {
    # If  mu  is a scalar, force it to be the same same length as  y
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  
  
  # Checking the length of  phi
  if ( length(phi) > 1) {
    # If  phi  not a scalar, check it is the same length as  y
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as the first input.\n")
    inputs_Error <- TRUE
    
  } else {
    # If  phi  is a scalar, force it to be the same same length as  y
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  # NOTE: The length of xi/power should have been checked in sort_Notation(), 
  #        so does not need checking here again.

  return( list(mu = mu, 
               phi = phi,
               inputs_Error = inputs_Error) )
} 




################################################################################

special_Cases <- function(y, mu, phi, power){
  # Special cass may be one of two types:
  # - based on the value of p:
  #   - p = 0: use Normal distribution
  #   - p = 1: use Poisson distribution
  #   - p = 2: use gamma distribution
  #   In this case, special_p_Cases is a scalar and is TRUE
  #
  # - other values of p, and hence based on value of y:
  #   - y < 0
  #   - y == 0
  #   In this case, special_y_Cases is a vector, and is TRUE when appropriate
  
  
  f <- array( 0, 
              dim = length(y) )
  special_p_Cases <- FALSE        # TRUE if special cases are defined by special values of p
  special_y_Cases <- rep(FALSE,   # TRUE where special cases are defined by special values of y
                         length(y) )
  
  # Special cases BASED ON VALUE OF p
  if ( (power == 0 ) | (power == 1) | (power== 2) ){
    # Special cases based on the value of p  
    special_p_Cases = TRUE
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
  } else {
    # Special cases BASED ON THE VALUES OF y (when y <= 0)
    special_y_Cases <- (y <= 0)
    if (any(special_y_Cases)) {

      # NEGATIVE VALUES
      y_Negative <- (y < 0)
      if (any(y_Negative) ) f[y_Negative] <- 0

      # ZERO VALUES
      y_Zero <- (y == 0)
      if (any(y_Zero)) {
        if ( (power > 0) & (power < 2) ) {
          f[y_Zero] <- exp( find_Lambda(mu[y_Zero], phi[y_Zero], power) )
        } else {
          f[y_Zero] <- 0  
        }
      }
    }
  }

  return( list(f = f,                                # vector
               special_p_Cases = special_p_Cases,    # scalar
               special_y_Cases = special_y_Cases) )  # vector
}


################################################################################

find_Lambda <- function(mu, phi, p){
  # Computes the value of lambda, such that P(Y = 0 ) = exp( -lambda) when 1 < p < 2
  if (p >= 2) stop("ERROR: lambda only sensible for  1 < xi < 2.")
  mu ^ (2 - p) / (phi * (2 - p) )
  
}