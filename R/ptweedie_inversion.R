#' @title Fourier Inversion Evaluation for the Tweedie Distribution Function
#' @name ptweedie_inversion
#' @description
#' Evaluates the distribution function (\acronym{df}) for Tweedie distributions 
#' using Fourier inversion, for given values of the dependent variable 
#' \code{y}, the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be in the case of evaluation problems.
#'
#' @usage ptweedie_inversion(q, mu, phi, power, lower.tail = TRUE, log.p = FALSE, IGexact = TRUE, 
#'                           verbose = FALSE, details = FALSE)
#'
#' @param q vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter.
#' @param phi the dispersion parameter.
#' @param lower.tail logical; if \code{TRUE} (the default) computes the distribution function \eqn{F(y)}; if \code{FALSE}, computes \eqn{1 - F(y)}. 
#' @param log.p logical; if \code{TRUE}, probabilities are returned as \eqn{\log(p)}, computed
#'   directly from the already-correctly-computed tail probability rather than as
#'   \code{log(ptweedie_inversion(...))} after the fact. The default is \code{FALSE}.
#' @param IGexact logical; if \code{TRUE} (the default), evaluate the inverse Gaussian distribution using the 'exact' values, otherwise uses inversion.
#' @param verbose logical; if \code{TRUE}, displays some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution and some information about the integration. The default is \code{FALSE}.
#' 
#' @return If \code{details = FALSE}, a numeric vector of the distribution function values; if \code{details = TRUE}, a list containing \code{CDF} (a vector of the values of the distribution function), \code{regions} (a vector of the number of integration regions used), and \code{exitstatus} (a vector, where a \code{1} for any value means a computational problem or target relative accuracy not reached, for the corresponding observation).
#' 
#' For special cases of \eqn{p} (i.e., \eqn{p = 0, 1, 2, 3}), where no inversion is needed, \code{regions} is set to \code{NA} for all values of \code{q}.
#' For special cases of \code{q} for other values of \eqn{p} (i.e., \eqn{P(Y = 0)}), \code{regions} is set to \code{NA}.
#'
#' @note
#' The 'exact' values for the inverse Gaussian distribution are not really exact, 
#' but evaluated using statsmod::[pdpqr]invgauss, all of which are very accurate
#' (Giner & Smyth, 2016).
#' 
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#' 
#' Giner G., Smyth G. K. (2016). 
#' statmod: probability calculations for the inverse Gaussian distribution. 
#' \emph{The R Journal},
#' \bold{8}(1), 339--351. 
#' \doi{doi:10.32614/RJ-2016-024}
#'
#' @examples
#' # Plot a Tweedie distribution function
#' y <- seq(0.01, 4, length = 50)
#' Fy <- ptweedie_inversion(y, mu = 1, phi = 1, power = 1.1)
#' plot(y, Fy, type = "l", lwd = 2, ylab = "Distribution function")
#' 
#' @keywords distribution
#' 
#' @export
ptweedie_inversion <- function(q, mu, phi, power, lower.tail = TRUE, log.p = FALSE, IGexact = TRUE,
                               verbose = FALSE, details = FALSE){ 
  ### NOTE: No notation checks
  
  # Check
  if (length(q) == 0L) {
    return(numeric(0))
  }
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(q, mu, phi, power)
  mu  <- out$mu
  phi <- out$phi

  # cdf    is the whole vector; the same length as  q.
  # All is resolved in the end.
  cdf <- numeric(length = length(q) )
  regions <- rep(NA, length(q)) 
  exitstatus_out <- integer(length(q))
  # exitstatus: 0 for values computed exactly (special cases); filled from the
  # Fortran for the rest. Always the same length as  q.
  
  # IDENTIFY SPECIAL CASES
  special_y_cases <- rep(FALSE, 
                         length(q) )
  if (verbose) cat("- Checking for special cases\n")
  out <- special_cases(q, mu, phi, power,
                       IGexact    = IGexact,
                       type       = "CDF",
                       lower.tail = lower.tail,
                       log.p      = log.p)
  
  
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases
  
  if (verbose & special_p_cases) cat("  - Special case for p used\n")
  if ( any(special_y_cases) ) {
    special_y_cases <- out$special_y_cases  
    if (verbose) cat("  - Special cases for first input found\n")
    cdf <- out$f # This is the final vector of results to return, filled with the special-case info.
    # NOTE: regions filled with zeros by default, so regions = 0 in these cases
  }
  
  if ( special_p_cases ) {
    cdf <- out$f
  } else {
    # NOT special p case; ONLY special y cases 
    
    # Now use FORTRAN on the remaining values:
    N_nonSpecial <- length(q) - sum(out$special_y_cases) 
      
    
    
    ### BEGIN: SET UP
    pSmall  <- ifelse( (power > 1) & (power < 2),
                       TRUE, 
                       FALSE )

    ### END SET UP
    
    ### Use re-scaling identity:
    ###  F(y; mu, phi) = F(y/mu; mu = 1, phi = mu^(p-2) * phi)
    mu_F  <- rep(1, length(q))
    q_F   <- q / mu
    phi_F <- mu^(power - 2) * phi
    
    # Scaling is numerically unsafe for extremely small q_F
    cut_off <- 1e-307
    use_scaled <- !special_y_cases & (q_F >= cut_off)
    use_direct <- !special_y_cases & (q_F <  cut_off)
    
    
    # CALL FORTRAN ROUTINES
    if (any(use_scaled)) {
      tmp <- .C(
        "twcomputation",
        N          = as.integer(sum(use_scaled)),
        power      = as.double(power),
        phi        = as.double(phi_F[use_scaled]),
        y          = as.double(q_F[use_scaled]),
        mu         = as.double(mu_F[use_scaled]),
        verbose    = as.integer(verbose),
        pdf        = as.integer(0),               # 0: FALSE, as this is the PDF
        tail       = as.integer(!lower.tail),     # 'tail' in the FORTRAN flags the UPPER tail (1 = Pr(Y > q)), hence !lower.tail
          # THE OUTPUTS:
        funvalue   = numeric(sum(use_scaled)),
        exitstatus = integer(sum(use_scaled)),
        relerr     = numeric(sum(use_scaled)),
        its        = integer(sum(use_scaled)),
        PACKAGE    = "tweedie"
      )
      
      cdf[use_scaled] <- tmp$funvalue
      regions[use_scaled] <- tmp$its
      exitstatus_out[use_scaled] <- tmp$exitstatus
    }
    if (any(use_direct)) {
      tmp <- .C(
        "twcomputation",
        N          = as.integer(sum(use_direct)),
        power      = as.double(power),
        phi        = as.double(phi[use_direct]),
        y          = as.double(q[use_direct]),
        mu         = as.double(mu[use_direct]),
        verbose    = as.integer(verbose),
        pdf        = as.integer(0),               # 0: FALSE, as this is the PDF
        tail       = as.integer(!lower.tail),     # 'tail' in the FORTRAN flags the UPPER tail (1 = Pr(Y > q)), hence !lower.tail
          # THE OUTPUTS:
        funvalue   = numeric(sum(use_direct)),
        exitstatus = integer(sum(use_direct)),
        relerr     = numeric(sum(use_direct)),
        its        = integer(sum(use_direct)),
        PACKAGE    = "tweedie"
      )
      
      cdf[use_direct] <- tmp$funvalue
      regions[use_direct] <- tmp$its
      exitstatus_out[use_direct] <- tmp$exitstatus
    }    

  }
  
  # CDF at this point is already the directly-computed
  # requested tail (per `tail`/lower.tail above), so log() here is safe and does not
  # involve any (1 - p) style subtraction.
  
  if ( special_p_cases ) {
    cdf <- out$f
      # cdf is already on the correct scale (log, if log.p was requested) --
      # special_cases() was called with log.p above, so do NOT log() again here.
  } else {
      # special_y_cases entries (if any) are already correctly scaled via
      # special_cases(); only the Fortran-computed entries are linear-scale
      # and need log() applied now.
    fortran_idx <- which(!special_y_cases)
    if (log.p && length(fortran_idx) > 0) {
      cdf[fortran_idx] <- log(cdf[fortran_idx])
    }
  }
  
  n_bad <- sum(exitstatus_out == 1L)
  if (n_bad > 0) {
    warning("ptweedie_inversion: the numerical integration did not reach the target accuracy for ", n_bad,
            " of ", length(exitstatus_out), " value(s); these may be inaccurate ",
            "(use details = TRUE to see which: exitstatus = 1).", call. = FALSE)
  }

  if (details) {
    return( list( cdf        = cdf,
                  regions    = regions,
                  exitstatus = exitstatus_out))
  } else {
    return(cdf)
  }
}

#' @rdname ptweedie_inversion
#' @export
ptweedie.inversion <- function(q, power, mu, phi, verbose, details){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "ptweedie.inversion()", 
                            with = "ptweedie_inversion()")
  ptweedie_inversion(q       = q, 
                     power   = power,
                     mu      = mu, 
                     phi     = phi, 
                     verbose = FALSE, 
                     details = FALSE)
}

