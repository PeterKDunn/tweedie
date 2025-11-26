#' Fourier Inversion Evaluation for the Tweedie Probability Function
#'
#' @description
#' Evaluates the probability density function (\acronym{pdf}) for Tweedie distributions using Fourier inversion, 
#' for given values of the dependent variable \code{y}, the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be used in the case of evaluation problems.
#'
#' @usage dtweedie_inversion(y, power, mu, phi, method = 3, verbose = FALSE, details = FALSE)
#'
#' @param y vector of quantiles.
#' @param power scalar; the power parameter \eqn{p}{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param method the method to use; one of \code{1}, \code{2}, or \code{3} (the default).
#' @param verbose logical; if \code{TRUE}, display some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, return a list with basic details of the integration. The default is \code{FALSE}.

#' @return A numeric vector of densities if \code{details=FALSE}; if \code{details=TRUE}, return a list with \code{density} (the density values) and \code{regions} (the number of integration regions used). 
#' 
#' @note
#' The three methods are described in Dunn & Smyth (2008).
#' 
#' @aliases dtweedie.inversion
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#' 
#' @export
dtweedie_inversion <- function(y, power, mu, phi, method = 3, verbose = FALSE, details = FALSE){ 
  ### NOTE: No checking of inputs
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(y, mu, phi, power)

  mu <- out$mu
  phi <- out$phi
  # density2 is the whole vector.
  # density is just the part where !special_y_cases.
  # All is resolved in the end
  density <- array(0,
                   dim = length(q) )
  if (details) regions <- array(0, dim = length(q))



  ### BEGIN SET UP
  N <- as.integer( length(y) )
  pSmall  <- ifelse( (power > 1) & (power < 2), 
                     TRUE, FALSE )
  # Initialise
  exitstatus_scalar <- as.integer(0)
  relerr_scalar     <- as.double(0.0)
  its_scalar        <- as.integer(0)
  ### END SET UP
  
  
  # Establish which method to use
  if ( is.null(method)){
    method <- array( dim = N)
  } else {
    method <- array( method, 
                     dim = N)
  }
  # There are three approaches ('method'), each a product of a simple bit
	# and a complicated bit computed in FORTRAN
  #
  # The methods are documented in Dunn and Smyth (2008):
  # - Method 1: Evaluate a(): compute a(y, phi) = f(y; 1, phi)
  # - Method 2: Rescale the mean to 1
  # - Method 3: Rescale y to 1 and evaluate b().
  #
  # If no method is explicitly requested, find the notional "optimal" method for each i.
  
 	### BEGIN ESTABLISH METHOD
  theta <- ( mu ^ (1 - power) - 1 ) / ( 1 - power )
  if ( ( abs(power - 2 ) ) < 1.0e-07 ){
    kappa <- log(mu) + (2 - power) * ( log(mu) ^ 2 ) / 2
  } else {
    kappa <- ( mu ^ (2 - power) - 1 ) / (2 - power)
  }

  # Method 1
  m1 <- exp( (y * theta - kappa ) / phi )
  dev <- tweedie_dev(y = y, 
                     mu = mu,
                     power = power )
  
  # Method 2
  m2 <- 1 / mu

  # Method 3
  m3 <- exp( -dev/(2 * phi) ) / y
    
  # Select method
  method_List <- array(c(m1, m2, m3), 
                       dim = c(length(y), 3))
  optimal_Method <- apply(method_List, 
                          MARGIN = 1, 
                          FUN = which.min)
  
  
  ### BEGIN: Set parameters for FORTRAN call, depending on method
  # mu = 1 for all methods:
  mu_F <- rep(1, length(y) )
  
  # Set up empty vector to fill for other methods:
  phi_F <- phi
  y_F <- y

  # Method 1 just uses the given  y  and  phi
  if (any(optimal_Method == 2)){
    use_M2 <- optimal_Method==2
    phi_F[ use_M2 ] <- phi[use_M2] / mu[use_M2] ^ (2 - power)
    y_F[ use_M2 ] <- y[use_M2]/mu[use_M2]
  }
  if (any(optimal_Method == 3)){
    use_M3 <- optimal_Method==3
    phi_F[ use_M3 ] <- phi[use_M3] / y[use_M3] ^ (2 - power)
    y_F[ use_M3 ] <- 1
  }
  ### END: Set parameters for FORTRAN call, depending on method

  tmp <- .C("twcomputation",
            N          = as.integer(N),
            power      = as.double(power),
            phi        = as.double(phi_F),
            y          = as.double(y_F),
            mu         = as.double(mu_F),
            verbose    = as.integer( verbose ),
            pdf        = as.integer(1),          # 1: TRUE, as this is the PDF
            # THE OUTPUTS:
            funvalue   = as.double(rep(0, N)),   # funvalue
            exitstatus = as.integer(0),          # exitstatus
            relerr     = as.double(0),           # relerr
            its        = as.integer(rep(0, N)),  # its
            PACKAGE    = "tweedie")
  
  den <- tmp$funvalue

  # Reconstruct
  if (any(optimal_Method == 1)){
    use_M1 <- optimal_Method==1
    density[use_M1] <- den[use_M1] * m1[use_M1]
  }    
  if (any(optimal_Method == 2)){
    density[use_M2] <- den[use_M2] * m2[use_M2]
  }  
  if (any(optimal_Method == 3)){
    density[use_M3] <- den[use_M3] * m3[use_M3]
  }

  # Return
  if (details) {
    return( list( density = density,
                  regions = tmp$its,
                  method = optimal_Method))
  } else {
    return(density)
  }
}


#' @export
dtweedie.inversion <- function(y, power, mu, phi, method = 3, verbose, details){ 
  .Deprecated("dtweedie_inversion", package = "tweedie")
  dtweedie_inversion(y = y, 
                     power = power, 
                     mu = mu, 
                     phi = phi, 
                     method = method, 
                     verbose = FALSE, 
                     details = FALSE)
}

