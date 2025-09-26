sort_Notation <- function(xi = NULL, power = NULL){
  
  if ( length(xi) > 1) stop(" The Tweedie index paramerer (xi) must be a single value.")
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