
#' @noRd
sort_notation <- function(xi = NULL, power = NULL){
  # Sorts out whether the  xi  or the  p/power  notation is being used,
  # so output presented appropriately
  
  if ( length(xi) > 1) stop(" The Tweedie index parameter (xi) must be a single value.")
  if ( length(power) > 1) stop(" The Tweedie index parameter (power) must be a single value.")
  
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

#' @noRd
check_inputs <- function(y, mu, phi, power, type = "standard"){
  # Checks that the inputs satisfy the necessary criteria (e.g., mu > 0).
  # Ensures that y, mu and phi are all vectors of the same length.

  # Since rtweedie(n, ...) does not have a first input as a vector, 
  # unlike [dpq]tweedie, flag this case and treat separately
  vector_Length <- ifelse(type == "random", 
                          max( length(mu), length(phi) ),
                          length(y))
  
  inputs_Error <- FALSE
  ### CHECKING VALUES ARE OK
  # Checking the input values for qtweedie
  if ( type == "quantile" ) {
    if ( any(y < 0) | any(y > 1) ) {
      stop("quantiles must be between 0 and 1.\n")
      inputs_Error <- TRUE
    }
  }
  
  ### CHECKING VALUES ARE OK
  # Checking the input values: n for rtweedie
  if ( type == "random" ) {
    y <- round(y)
    if ( y <= 0 ) {
      stop("n must be a positive integer.\n")
      inputs_Error <- TRUE
    }
  }

  # Checking the input values: power
  if ( any(power < 1 ) ) {
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


  
  
  ### CHECKING LENGTHS ARE OK
  if (type == "random") {
    mu  <- rep_len(mu, y)
    phi <- rep_len(phi, y)
  } else {
    # Checking the length of  mu
    if ( length(mu) > 1) {
      # If  mu  not a scalar, check it is the same length as  y
      if ( length(mu) != length(y) ) {
        stop("mu must be scalar, or the same length as the first input.\n")
        inputs_Error <- TRUE
      }
    } else {
      # If  mu  is a scalar, force it to be the same same length as  y
      mu <- rep_len(mu, length(y) )
      # A vector of all mu's
    }
    
    
    # Checking the length of  phi
    if ( length(phi) > 1) {
      # If  phi  not a scalar, check it is the same length as  y
      if ( length(phi) != length(y) ) {
        stop("phi must be scalar, or the same length as the first input.\n")
        inputs_Error <- TRUE
      }
    } else {
      # If  phi  is a scalar, force it to be the same same length as  y
      phi <- rep_len(phi, length(y) )
      # A vector of all phi's
    }
  }
  
  # NOTE: The length of xi/power should have been checked in sort_notation(), 
  #        so does not need checking here again.

  return( list(mu = mu, 
               phi = phi,
               inputs_Error = inputs_Error) )
} 




################################################################################

#' @noRd
special_cases <- function(y, mu, phi, power, type = "PDF", verbose = FALSE, IGexact = TRUE, 
                          lower.tail = TRUE, log.p = FALSE){
  # Special cases may be one of two types:
  # - based on the value of p:
  #   - p = 0: use Normal distribution
  #   - p = 1: use Poisson distribution
  #   - p = 2: use gamma distribution
  #   - p = 3: use inverse Gaussian distribution, via
  #   In this case, special_p_cases is a scalar and is TRUE
  #
  # - other values of p, and hence based on value of y:
  #   - y < 0
  #   - y == 0
  #   In this case, special_y_cases is a vector, and is TRUE when appropriate
  

  f <- numeric( length(y) )
  special_p_cases <- FALSE        # TRUE if special cases are defined by special values of p: SCALAR
  special_y_cases <- rep(FALSE,   # TRUE where special cases are defined by special values of y: VECTOR
                         length(y) )
  
  # Special cases BASED ON VALUE OF p
  if ( (power == 0 ) | (power == 1) | (power == 2) | (power == 3)){
    if (verbose) cat("Special cases in p found ")
    # Special cases based on the value of p  
    
    special_p_cases <- TRUE
    
    # CASE: Normal (p=0)
    if ( power == 0) {
      if (verbose) cat("power = 0 (normal case)\n")
      if (type == "PDF") {
        f <- stats::dnorm( y, 
                           mean = mu, 
                           sd   = sqrt(phi))
      } else {
        f <- stats::pnorm( y, 
                           mean       = mu, 
                           sd         = sqrt(phi), 
                           lower.tail = lower.tail,
                           log.p      = log.p)
      }
    }
    
    # CASE: Poisson (p=1)
    if ( power == 1) {
      if (verbose) cat("power = 1 (Poisson case)\n")
      if (type == "PDF"){
        f <- stats::dpois(y/phi, 
                          lambda = mu / phi )
      } else {
        f <- stats::ppois(y/phi, 
                          lambda     = mu / phi,
                          lower.tail = lower.tail,
                          log.p      = log.p)
      }
    }
    
    # CASE: gamma (p=2)
    if ( power == 2 ) {
      if (verbose) cat("power = 2 (gamma case)\n")
      if (type == "PDF") {
        f <- stats::dgamma( y,
                            scale = mu * phi, 
                            shape = 1 / phi)
      } else {
        f <- stats::pgamma( y,
                            scale      = mu * phi, 
                            shape      = 1 / phi, 
                            lower.tail = lower.tail,
                            log.p      = log.p)
      }
    }
    
    # CASE: inverse Gaussian (p=3)
    if (power == 3) {
      if (IGexact) {
        if (verbose) cat("power = 3 (inverse Gaussian case)\n")
        if (type == "PDF") {
          f <- statmod::dinvgauss(x          = y, 
                                  mean       = mu, 
                                  dispersion = phi)
        } else {
          f <- statmod::pinvgauss(q          = y, 
                                  mean       = mu, 
                                  dispersion = phi, 
                                  lower.tail = lower.tail,
                                  log.p      = log.p)
        }
      } else {
        special_p_cases = FALSE
      }
    }

  } else {
    # Special cases BASED ON THE VALUES OF y (when y <= 0)
    if (verbose) cat("Special cases in y found.\n")
    
    special_y_cases <- (y <= 0)
    if (any(special_y_cases)) {
      # NEGATIVE VALUES
      y_Negative <- (y < 0)
      # P(Y <= y) = 0 and P(Y > y) = 1 for y < 0
      if (any(y_Negative) ) f[y_Negative] <- if (lower.tail) { if (log.p) -Inf else 0 } else { if (log.p) 0 else 1 }
      y_Zero <- (y == 0)
      if (any(y_Zero)) {
        if ( (power > 0) & (power < 2) ) {
          lambda_z <- tweedie_lambda(mu[y_Zero], phi[y_Zero], power)
          if (lower.tail) {
            f[y_Zero] <- if (log.p) -lambda_z else exp(-lambda_z)
          } else {
            # S(0) = 1 - exp(-lambda); use log1p(-exp(-lambda)) directly for log.p,
            # not log(1 - exp(-lambda)) after forming the linear value
            f[y_Zero] <- if (log.p) {
              log1p(-exp(-lambda_z))
            } else {
              -expm1(-lambda_z)   # = 1 - exp(-lambda), computed via expm1 for precision when lambda_z is small
            }
          }
        } else {
          # p >= 2: no point mass at zero, so P(Y <= 0) = 0 and P(Y > 0) = 1
          f[y_Zero] <- if (lower.tail) { if (log.p) -Inf else 0 } else { if (log.p) 0 else 1 }
        }
      }
    }
  }

  return( list(f = f,                                # vector
               special_p_cases = special_p_cases,    # scalar
               special_y_cases = special_y_cases) )  # vector
}



################################################################################

#' @noRd
write_f90_gauss_table <- function(n, 
                                  digits           = 17,
                                  array_name_nodes = "absc",
                                  array_name_wts   = "wts",
                                  values_per_line  = 3,
                                  file             = NULL) {
  
  # USAGE:  write_f90_gauss_table(200, array_name_nodes = "absc200", array_name_wts = "wts200")
  # Written so I can easily explore using different numbers of quadrature nodes  
  
  gq <- statmod::gauss.quad(n, kind = "legendre")
  
  pos          <- gq$nodes > 0
  ord          <- order(gq$nodes[pos])
  nodes_half   <- gq$nodes[pos][ord]
  weights_half <- gq$weights[pos][ord]
  
  half_n <- length(nodes_half)
  
  fmt_f90_block <- function(x, array_name, half_n) {
    lits <- sprintf(paste0("%.", digits - 1, "fD00"), x)
    
    chunks <- split(lits, 
                    ceiling(seq_along(lits) / values_per_line))
    n_chunks <- length(chunks)
    
    body_lines <- character(n_chunks)
    for (i in seq_len(n_chunks)) {
      vals <- chunks[[i]]
      if (i < n_chunks) {
        # more values follow: comma after last value on this line, then continue
        body_lines[i] <- paste0("    ", paste(vals, collapse = ", "), ", &")
      } else {
        # last chunk: no comma after the final value, but STILL continue
        # onto the next line, since '/)' is on its own separate line
        body_lines[i] <- paste0("    ", paste(vals, collapse = ", "), " &")
      }
    }
    
    header <- sprintf("  REAL(KIND=C_DOUBLE), PARAMETER, DIMENSION(%d) :: %s = (/ &",
                      half_n, array_name)
    footer <- "  /)"
    c(header, body_lines, footer)
  }
  
  out <- c(
    fmt_f90_block(nodes_half,   
                  array_name_nodes, 
                  half_n),
    "",
    fmt_f90_block(weights_half, 
                  array_name_wts,   
                  half_n)
  )
  
  cat(paste(out, collapse = "\n"), "\n")
  
  if (!is.null(file)) {
    writeLines(out, file)
    message("Written to: ", file)
  }
  
  invisible(out)
}

################################################################################
