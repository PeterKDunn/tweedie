#' @title Tweedie distributions
#' @name Tweedie
#' @aliases dtweedie ptweedie qtweedie rtweedie
#'   
#' @description Density, distribution function, quantile function and random 
#' generation for the the Tweedie family of distributions, with mean \code{mu}, 
#' dispersion parameter \code{phi} and variance power \code{power} (or \code{xi}, 
#' a synonym for \code{power}).
#' 
#' @usage 
#' dtweedie(y, xi = NULL, mu, phi, power = NULL, verbose = FALSE, log = FALSE)
#' ptweedie(q, xi = NULL, mu, phi, power = NULL, lower.tail = TRUE, 
#'   log.p = FALSE, verbose = FALSE)
#' qtweedie(p, xi = NULL, mu, phi, power = NULL)
#' rtweedie(n, xi = NULL, mu, phi, power = NULL)
#' 
#' @details
#' Tweedie \acronym{edm}s belong to the class of exponential dispersion models 
#' (\acronym{edm}s), known for their role in generalized linear models (\acronym{glm}s). 
#' The Tweedie distributions are the \acronym{edm}s with a variance of the form 
#' \eqn{\mbox{var}[Y] = \phi\mu^p}{var[Y] = phi*mu^p} where \eqn{p \ge 1}{p >= 1}.
#' \emph{This function only evaluates for \eqn{p \ge 1}{p >= 1}.}
#'
#' Special cases are the Poisson (\eqn{p = 1} with \eqn{\phi = 1}{phi = 1}), 
#' gamma (\eqn{p = 2}), and inverse Gaussian (\eqn{p = 3}) distributions.
#' Evaluation is difficult for \eqn{p}{p} outside of \eqn{p = 0, 1, 2, 3}{power = 0, 1, 2, 3}. 
#' For the density, \code{dtweedie()} uses one of three methods for each value, 
#' depending on \eqn{p}{p} and \eqn{z = \phi y^{p - 2}}{z = phi * y^(p - 2)} 
#' (the dispersion after rescaling \eqn{y}{y} to one):
#' \enumerate{
#'   \item the saddlepoint approximation with its asymptotic correction, 
#'         where \eqn{z}{z} is small enough for its estimated relative error to be below \eqn{10^{-13}}{1e-13};
#'   \item otherwise, the series (\code{dtweedie_series}) for \eqn{1 < p < 2}{1 < power < 2}, and 
#'         for \eqn{p > 2}{power > 2} where its estimated relative error is below \eqn{10^{-12}}{1e-12};
#'   \item otherwise, the Fourier inversion (\code{dtweedie_inversion}).
#' }
#' Against high-precision reference values, the relative error of the density is 
#' below about \eqn{4 \times 10^{-13}}{4e-13}.
#' (Before version 3.2.0, \code{dtweedie()} also interpolated from stored values 
#' of the Fourier inversion, which was less accurate.)
#'
#' When \eqn{1<p<2}{1 < power < 2}, the density function include a positive 
#' probability for \eqn{Y = 0}.
#'
#' @section Note:
#' \code{dtweedie()} and \code{ptweedie()} are the only functions generally to be called by users. 
#' Consequently, all checks on the function inputs are only performed in these functions.
#'
#' @param y vector of quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken 
#'    to be the number required (as for \code{\link[stats]{rnorm}}).
#' @param xi scalar; the value of \eqn{\xi}{xi} such that the variance is 
#'    \eqn{\mbox{var}[Y]=\phi\mu^{\xi}}{var[Y] = phi * mu^xi}. A synonym for \code{power}.
#' @param mu vector of mean \eqn{\mu}{mu}.
#' @param phi vector of dispersion parameters \eqn{\phi}{phi}.
#' @param power scalar; a synonym for \eqn{\xi}{xi}, the Tweedie index parameter.
#' @param lower.tail logical; if \code{TRUE} (the default), computes the lower 
#'    tail \eqn{F(y)}; otherwise computes \eqn{S(y) = 1 - F(y)}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned as \eqn{\log(p)}.
#'    The default is \code{FALSE}.
#' @param verbose logical; if \code{TRUE}, some details of the algorithms used is
#'    shown. The default is \code{FALSE}.
#' @param log logical; if \code{TRUE}, \code{dtweedie()} returns the log of the density,
#'    computed without underflow (so it stays finite far into the tails, where the
#'    density itself underflows to zero). The default is \code{FALSE}.
#'
#' @return
#' \code{dtweedie()} gives the density, \code{ptweedie()} gives the distribution function,
#' \code{qtweedie()} gives the quantile function, 
#' and \code{rtweedie} generates random deviates.
#' 
#' As for R's own functions for distributions (such as \code{\link[stats]{dnorm}}), 
#' the result of \code{dtweedie()}, \code{ptweedie()} and \code{qtweedie()} has 
#' the same length and form as the first input (\code{y}, \code{q} or \code{p}): 
#' its attributes are kept, so (for example) a matrix gives a matrix, and a 
#' named vector gives a named vector.
#' (\code{mu} and \code{phi} must be scalars, or the same length as the first input.)
#' \code{rtweedie()} returns a vector of length \code{n}, without attributes.
#' 
#' @importFrom stats dgamma dnorm dpois
#' @importFrom graphics lines legend plot
#' 
#' @seealso \code{\link{dtweedie_series}}, \code{\link{dtweedie_inversion}}, 
#' \code{\link{ptweedie_series}}, \code{\link{ptweedie_inversion}}, \code{\link{dtweedie_saddle}}, \code{\link{tweedie_lambda}}
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#' 
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' Jorgensen, B. (1997).
#' \emph{Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' @examples
#' # Compute a Tweedie density
#' power <- 1.1
#' mu <- 1
#' phi <- 1
#' y <- seq(0, 5, by = 0.5)
#' dtweedie(y, power = power, mu = mu, phi = phi)
#' 
#' # Compare to the saddlepoint density
#' dtweedie_saddle(y = y, power = power, mu = mu, phi = phi)
#' 
#' # The DF:
#' ptweedie(y, power = power, mu = mu, phi = phi)
#' 
#' @keywords distribution
#' @export
dtweedie <- function(y, xi = NULL, mu, phi, power = NULL, verbose = FALSE, log = FALSE){
  # THE REGIMES (since 2026; no interpolation)
  #
  # Special values of p (p = 1, 2, 3; also p = 0) use the closed forms, and
  # y <= 0 is exact (the point mass Pr(Y = 0) = exp(-lambda) when 1 < p < 2).
  # Otherwise, for each y > 0, with the dispersion after rescaling y to 1,
  #
  #   z = phi y^(p - 2)
  #
  # (by the scaling property, the relative accuracy of each method depends on
  # p and z only), three methods are used:
  #
  #   1. The saddlepoint approximation with its asymptotic correction
  #      (dtweedie_saddle_corrected()), wherever its estimated relative error is
  #      below 1e-13: that is, wherever z is small enough. This needs z below
  #      about 0.07 for p from 1.1 to 2.05, 0.03 for p = 2.5, 3e-3 for p = 4,
  #      4e-4 for p = 8 and 1e-4 for p = 12 (for 1 < p < 2 the density also has
  #      ripples that the approximation misses, which matter for p near 1: there
  #      z must be below 7e-3 when p = 1.01). Its actual relative error there
  #      is about 1e-13 or smaller. (The uncorrected saddlepoint approximation
  #      has relative error of about z p |p - 3| / 24.)
  #
  #   2. Otherwise, for 1 < p < 2: the series (dtweedie_series()), a sum of
  #      positive terms, accurate to about 1e-13 wherever it is used here (z
  #      not small means few terms: about 1/{z (2 - p)} at the peak).
  #
  #   3. Otherwise, for p > 2: the series, if its estimated relative error is
  #      below 1e-12 (its terms alternate in sign, and it cancels badly once
  #      its peak index k_max = 1/{z (p - 2)} is more than a few times p - 1);
  #      otherwise the Fourier inversion (dtweedie_inversion(), method 3:
  #      y rescaled to 1), accurate to about 1e-12 or better there.
  #
  # If the inversion fails (it reports exitstatus = 1), whichever of the
  # series and the saddlepoint approximation has the smaller estimated error
  # is used instead, with a warning.
  #
  # Against high-precision reference values (p from 1.01 to 12, z from 1e-12
  # to 1e3, the body and both tails down to tail probabilities of 1e-30), the
  # largest relative error is about 4e-13.
  #
  # Before 2026, dtweedie() also interpolated from stored values of the
  # inversion, which was faster but much less accurate in places (relative
  # errors of about 1e-9 for most p, and about 2.5e-7 for 7 <= p <= 10).
  
  # (The thresholds are in tweedie_regime_settings(); dtweedie_regime() says
  # which regime each value is in.)
  S <- tweedie_regime_settings()
  
  ### BEGIN: preliminary work
  
  # SORT OUT THE NOTATION (i.e., xi VS power)
  if (verbose) cat("- Checking notation\n")
  out   <- sort_notation(xi = xi, power = power)
  xi    <- out$xi
  power <- out$power
  
  xi.notation    <- out$xi.notation
  index.par      <- out$index.par
  index.par.long <- out$index.par.long 
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(y, mu, phi, power)
  mu  <- out$mu
  phi <- out$phi
  
  # IDENTIFY SPECIAL CASES
  if (verbose) cat("- Checking for special cases\n")
  # (The log of the density is kept alongside the density throughout, from
  # each method's own log scale, so that log = TRUE does not underflow.)
  out <- special_cases(y, mu, phi, power,
                       type  = "PDF",
                       log.p = log)
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases
  if (log) {
    logdens <- out$f
    density <- exp(logdens)
  } else {
    density <- out$f
    logdens <- base::log(density)
  }
  if (verbose & special_p_cases) cat("  - Special case for p used\n")
  if (verbose & any(special_y_cases)) cat("  - Special cases for first input found\n")
  
  ### END: preliminary work
  
  if ( !special_p_cases ) {
    
    # The values still to compute: y > 0 (y = Inf has density zero)
    todo <- which( !special_y_cases & !is.na(y) )
    y_inf <- todo[ y[todo] == Inf ]
    density[y_inf] <- 0
    logdens[y_inf] <- -Inf
    todo <- setdiff(todo, y_inf)
    
    method <- rep(NA_character_, length(y))  # the method used for each value (for inspection; not returned)
    
    if ( length(todo) > 0 ) {
      
      regime <- dtweedie_regime(y     = y[todo], 
                                mu    = mu[todo], 
                                phi   = phi[todo], 
                                power = power)
      
      # 1. The corrected saddlepoint approximation, where accurate (a log
      #    density of -Inf is fine: the density underflows to zero)
      use_sp <- (regime == "saddlepoint")
      if ( any(use_sp) ) {
        sp <- dtweedie_saddle_corrected(y     = y[ todo[use_sp] ], 
                                        mu    = mu[ todo[use_sp] ], 
                                        phi   = phi[ todo[use_sp] ], 
                                        power = power)
        density[ todo[use_sp] ] <- sp$density
        logdens[ todo[use_sp] ] <- sp$logf
        method[ todo[use_sp] ]  <- "saddlepoint"
      }
      rest <- todo[!use_sp]
      if (verbose) cat("- Corrected saddlepoint approximation:", sum(use_sp), "value(s)\n")
      
      if ( length(rest) > 0 ) {
        
        if (power < 2) {
          # 2. 1 < p < 2: the series
          ser <- dtweedie_series(y       = y[rest], 
                                 mu      = mu[rest], 
                                 phi     = phi[rest],
                                 power   = power,
                                 details = TRUE)
          density[rest] <- ser$density
          logdens[rest] <- ser$logf
          method[rest]  <- "series"
          if (verbose) cat("- Series:", length(rest), "value(s)\n")
          
        } else {
          # 3. p > 2: the series where its estimated error is small...
          ser <- suppressWarnings( 
                   dtweedie_series(y       = y[rest], 
                                   mu      = mu[rest], 
                                   phi     = phi[rest],
                                   power   = power,
                                   details = TRUE) )
          use_ser <- !is.na(ser$relerr) & (ser$relerr < S$d_tol_series)
          density[ rest[use_ser] ] <- ser$density[use_ser]
          logdens[ rest[use_ser] ] <- ser$logf[use_ser]
          method[ rest[use_ser] ]  <- "series"
          if (verbose) cat("- Series:", sum(use_ser), "value(s)\n")
          
          # ... otherwise the inversion
          inv <- rest[!use_ser]
          if ( length(inv) > 0 ) {
            if (verbose) cat("- Inversion:", length(inv), "value(s)\n")
            inv_out <- suppressWarnings(
                         dtweedie_inversion(y       = y[inv],
                                            mu      = mu[inv],
                                            phi     = phi[inv],
                                            power   = power,
                                            method  = 3,
                                            details = TRUE) )
            density[inv] <- inv_out$density
            logdens[inv] <- inv_out$logdensity
            method[inv]  <- "inversion"
            
            # If the inversion failed: the series or the saddlepoint
            # approximation, whichever has the smaller estimated error
            failed <- (inv_out$exitstatus == 1L) | is.na(inv_out$density)
            if ( any(failed) ) {
              k      <- which(failed)
              re_ser <- ser$relerr[!use_ser][k]
              re_ser[is.na(re_ser)] <- Inf
              sp_k   <- dtweedie_saddle_corrected(y     = y[ inv[k] ], 
                                                  mu    = mu[ inv[k] ], 
                                                  phi   = phi[ inv[k] ], 
                                                  power = power)
              dens_sp_k <- sp_k$density
              re_sp  <- sp_k$relerr
              re_sp[is.na(re_sp) | !is.finite(dens_sp_k)] <- Inf
              pick_ser <- (re_ser <= re_sp) & is.finite(re_ser)
              pick_sp  <- !pick_ser & is.finite(re_sp)
              density[ inv[k][pick_ser] ] <- ser$density[!use_ser][k][pick_ser]
              logdens[ inv[k][pick_ser] ] <- ser$logf[!use_ser][k][pick_ser]
              method[ inv[k][pick_ser] ]  <- "series"
              density[ inv[k][pick_sp] ]  <- dens_sp_k[pick_sp]
              logdens[ inv[k][pick_sp] ]  <- sp_k$logf[pick_sp]
              method[ inv[k][pick_sp] ]   <- "saddlepoint"
              # (No warning where the density underflows to zero, whatever 
              # the method's relative error)
              n_warn <- sum( !is.na(density[ inv[k] ]) & (density[ inv[k] ] > 0) )
              if (n_warn > 0) {
                warning("dtweedie: the Fourier inversion did not reach the target accuracy for ", 
                        n_warn, " value(s); ",
                        "the series or the saddlepoint approximation was used for these instead, ",
                        "and they may be inaccurate.", 
                        call. = FALSE)
              }
            }
          }
        }
      }
    }
  }
  
  # Restoration of whole vector, and sanity fixes
  density[is.na(y)] <- NA
  logdens[is.na(y)] <- NA
  neg <- !is.na(density) & (density < 0)
  if (any(neg)) {
    density[neg] <- 0
    logdens[neg] <- -Inf
  }
  logdens[ is.nan(logdens) & !is.na(y) ] <- -Inf   # (the log of a negative value)
  density <- if (log) logdens else density

  # The result has the same form (dim, names, ...) as  y, as for dnorm() etc.
  return( match_shape(density, y) )

}

################################################################################

#' @export
ptweedie <- function(q, xi = NULL, mu, phi, power = NULL, lower.tail = TRUE, log.p = FALSE, 
                     verbose = FALSE){
  ### BEGIN preliminary work
  
  # SORT OUT THE NOTATION (i.e., xi VS power)
  if (verbose) cat("- Checking notation\n")
  out <- sort_notation(xi = xi, power = power)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long 
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(q, mu, phi, power)
  mu <- out$mu
  phi <- out$phi
  
  # MISSING VALUES: NA in, NA out; the other values are computed as usual
  # (before, an NA gave an error)
  if ( any(is.na(q)) ) {
    ok  <- !is.na(q)
    ans <- rep(NA_real_, length(q))
    if ( any(ok) ) {
      ans[ok] <- ptweedie(q          = q[ok], 
                          mu         = mu[ok], 
                          phi        = phi[ok], 
                          power      = power, 
                          lower.tail = lower.tail, 
                          log.p      = log.p, 
                          verbose    = verbose)
    }
    return( match_shape(ans, q) )
  }
  
  f <- array(0,
             dim = length(q) )
  
  # IDENTIFY SPECIAL CASES
  special_y_cases <- rep(FALSE, length(q))
  if (verbose) cat("- Checking for special cases\n")
  out <- special_cases(q, mu, phi, power, 
                       type       = "CDF",
                       lower.tail = lower.tail,
                       log.p      = log.p)
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases
  
  if (verbose & special_p_cases) cat("  - Special case for p used\n")
  if ( any(special_y_cases) ) {
    special_y_cases <- out$special_y_cases  
    if (verbose) cat("  - Special cases for first input found\n")
    f <- out$f
  }
  
  ### END preliminary work
  
  if ( special_p_cases ) {
    f <- out$f
  } else {
    # NOT special p case; ONLY special y cases 
    
    if ( power > 2 ) {
      # For p > 2, with xi = mu^(2 - p) / {phi (p - 2)} and the dispersion
      # after scaling to mean 1, phis = phi mu^(p - 2), three regimes:
      #
      #   xi <= 1e4        : the inversion; tail probabilities below 1e-5 are
      #                      recomputed by integrating the density (between
      #                      1e-10 and 1e-5, only if that agrees with the
      #                      inversion to 1e-12, since the inversion is then the
      #                      more reliable of the two).
      #   1e4 < xi <= 1e7  : the inversion's error grows with xi (relative
      #                      errors of about 1e-8 at tail probabilities of 1e-4
      #                      when xi = 1e5, and 3e-6 when xi = 1e7, for every p
      #                      tested), while integrating the density stays
      #                      accurate (1e-9 or better, against 240-bit reference
      #                      values, for p from 2.001 to 8, xi up to 1e7, and
      #                      tail probabilities from 1e-2 to 1e-33). So tails
      #                      below 1e-2 are integrated, and the integrated value
      #                      is used whenever the integration succeeds. In the
      #                      body, the inversion is used for xi <= 1e5; for
      #                      1e5 < xi, the density is integrated in the body
      #                      too (the inversion's error in the body reaches
      #                      about 1e-8 near p = 2 at xi = 1e7, against about
      #                      1e-9 for the integration; the integration is
      #                      slower, about 0.1 s per value), with the inversion
      #                      where the integration returns no value. For
      #                      1e6 < xi and phis <= 1e-5, the saddlepoint
      #                      approximation is used in the body instead (relative
      #                      error below about 3e-10).
      #   xi > 1e7         : the saddlepoint approximation, for every y; and
      #                      also at xi = 1e7 when phis <= 1e-5 (then p >= 2.01),
      #                      where it is more accurate than the integration. Its
      #                      error falls as xi grows (about 2e-8 at xi = 1e7 for
      #                      p = 2.001, and below 1e-9 there for p >= 2.01),
      #                      while the integration can fail outright by
      #                      xi = 1e8.
      #
      # Where the inversion fails (for xi <= 1e4, at very small y when p is near
      # 2, it can return about 0.5) and none of the above replaces its value,
      # the density is integrated; failing that, the saddlepoint approximation
      # is used if it gives a probability (rough there, but far better than the
      # failed value). Such values are still flagged by a warning.
      #
      # The saddlepoint approximation's accuracy depends on phis, not on xi:
      # near p = 2, xi is large even when phis is not small (xi = 1e5 when
      # p = 2.001 and phis = 0.01), and there it is much less accurate (errors
      # of 1e-6 to 1e-4), which is why it is not simply used for large xi.
      # Where it returns no value (NA), the other methods are used.
      # (The thresholds are in tweedie_regime_settings(); see also
      # ptweedie_regime(), which says which regime each value is in.)
      S       <- tweedie_regime_settings()
      use_inv <- !special_y_cases
      xi_d    <- mu^(2 - power) / (phi * (power - 2))
      phis    <- phi * mu^(power - 2)
      if ( any(use_inv) ) {
        # (xi within rounding of 1e7 counts as xi = 1e7: it is included only
        # when phis <= 1e-5; near p = 2, at xi = 1e7, the integration is the
        # more accurate)
        use_sp <- use_inv & ptweedie_use_saddle(xi_d, phis)
        if (any(use_sp)) {
          if (verbose) cat("- With p > 2: using the saddlepoint approximation for", sum(use_sp), "value(s) (very small dispersion)\n")
          f_sp <- ptweedie_saddle(q          = q[use_sp],
                                  mu         = mu[use_sp],
                                  phi        = phi[use_sp],
                                  power      = power,
                                  lower.tail = lower.tail,
                                  log.p      = log.p)
          f[use_sp] <- f_sp
          use_inv[which(use_sp)[!is.na(f_sp)]] <- FALSE
        }
      }
      if ( any(use_inv)) {
        if (verbose) cat("- With p > 2: use inversion\n")

        f_TMP <- suppressWarnings(
                   ptweedie_inversion(q          = q[use_inv],
                                      mu         = mu[use_inv],
                                      phi        = phi[use_inv],
                                      power      = power,
                                      lower.tail = lower.tail,
                                      log.p      = log.p,
                                      verbose    = verbose,
                                      details    = TRUE) )
        f[use_inv] <- f_TMP$cdf
        not_conv <- (f_TMP$exitstatus == 1L)
        
        # The inversion's error is ABSOLUTE, so a tail probability P has a
        # relative error of about (absolute error)/P. Small tails are
        # therefore recomputed by integrating the density, which keeps its
        # relative accuracy far into the tails (see ptweedie_tail_integrate;
        # Gauss-Laguerre quadrature first, so this is usually cheap). See the
        # regimes above for the thresholds.
        tail_switch <- S$p_tail_small              # for xi <= 1e4
        idx      <- which(use_inv)
        trust    <- xi_d[idx] > S$p_xi_trust       # 1e4 < xi <= 1e7
        tail_lin <- if (log.p) exp(f_TMP$cdf) else f_TMP$cdf
        small    <- is.na(tail_lin) | (tail_lin < ifelse(trust, S$p_tail_trust, tail_switch))
        not_rel  <- rep(FALSE, length(idx))
        if (any(small)) {
          if (verbose) cat("- Small tail probabilities:", sum(small), "value(s) recomputed by integrating the density\n")
          for (k in which(small)) {
            j <- idx[k]
            # For xi <= 1e4, between 1e-10 and 1e-5 only the cheap
            # Gauss-Laguerre rules are tried, since the inversion value is a
            # good fallback there
            very_small <- is.na(tail_lin[k]) || (tail_lin[k] < S$p_tail_tiny)
            ti <- ptweedie_tail_integrate(q          = q[j], 
                                          mu         = mu[j], 
                                          phi        = phi[j],
                                          power      = power, 
                                          lower.tail = lower.tail,
                                          fallback   = very_small || trust[k])
            if (ti$ok && !very_small && !trust[k] &&
                !(abs(exp(ti$logp) - tail_lin[k]) < S$p_agree)) {
              next   # disagrees with the inversion: keep the inversion value
            }
            if (ti$ok) {
              f[j] <- if (log.p) {
                ti$logp 
              } else {
                 exp(ti$logp)
              }
              not_conv[k] <- FALSE
              # A tail so small that the density underflows: the integration
              # returns log P = -Inf. For 1e4 < xi the saddlepoint
              # approximation (relative error below 1e-4 there) gives the
              # finite log instead.
              if (log.p && trust[k] && (ti$logp == -Inf)) {
                ls <- ptweedie_saddle(q = q[j], mu = mu[j], phi = phi[j], power = power,
                                      lower.tail = lower.tail,
                                      log.p      = TRUE)
                if (!is.na(ls)) f[j] <- ls
              }
            } else if (very_small) {
              if (is.na(tail_lin[k])) {
                # No value at all (the inversion returned NA, and it could
                # not be recomputed): report it as not computed, not as
                # 'accurate only in absolute terms'
                not_conv[k] <- TRUE
              } else {
                not_rel[k] <- TRUE
              }
            }
          }
        }

        # The body, for 1e6 < xi <= 1e7 with phis <= 1e-5: the saddlepoint
        # approximation is more accurate than the inversion there
        body_sp <- !small & trust & (xi_d[idx] > S$p_xi_body_sp) & (phis[idx] <= S$p_phis_saddle)
        if (any(body_sp)) {
          j    <- idx[body_sp]
          f_sp <- ptweedie_saddle(q          = q[j],
                                  mu         = mu[j],
                                  phi        = phi[j],
                                  power      = power,
                                  lower.tail = lower.tail,
                                  log.p      = log.p)
          ok_sp <- !is.na(f_sp)
          f[j[ok_sp]] <- f_sp[ok_sp]
          not_conv[which(body_sp)[ok_sp]] <- FALSE
          body_sp[which(body_sp)[!ok_sp]] <- FALSE
        }

        # The body, for 1e5 < xi <= 1e7 (where the saddlepoint is not used):
        # integrating the density is more accurate than the inversion there.
        # Where the integration returns no value, the inversion's is kept.
        body_int <- !small & (xi_d[idx] > S$p_xi_body_int) & !body_sp
        if (any(body_int)) {
          if (verbose) cat("- Body values with small dispersion:", sum(body_int), "value(s) recomputed by integrating the density\n")
          for (k in which(body_int)) {
            j  <- idx[k]
            ti <- ptweedie_tail_integrate(q          = q[j],
                                          mu         = mu[j],
                                          phi        = phi[j],
                                          power      = power,
                                          lower.tail = lower.tail,
                                          fallback   = TRUE)
            if (ti$ok) {
              f[j] <- if (log.p) ti$logp else exp(ti$logp)
              not_conv[k] <- FALSE
            }
          }
        }

        # A flagged value close to one (e.g. P(Y > y) for y far below the
        # mean, when the inversion could not converge): its absolute error
        # is what matters, so compute instead the OTHER tail, which is then
        # small, by integrating the density, and return one minus it.
        near_one <- not_conv & !not_rel & !is.na(tail_lin) & (tail_lin > 1 - tail_switch)
        for (k in which(near_one)) {
          j  <- idx[k]
          ti <- ptweedie_tail_integrate(q          = q[j],
                                        mu         = mu[j],
                                        phi        = phi[j],
                                        power      = power,
                                        lower.tail = !lower.tail,
                                        fallback   = TRUE)
          if (ti$ok && (exp(ti$logp) < tail_switch)) {
            f[j] <- if (log.p) {
              log1p(-exp(ti$logp)) 
            } else {
              -expm1(ti$logp)
            }
            not_conv[k] <- FALSE
          }
        }

        # Where the inversion failed in what it took to be the body, and
        # nothing above replaced its value (for xi <= 1e4 this happens at very
        # small y when p is near 2, where the inversion can return about 0.5
        # whatever the true value): integrate the density; failing that, use
        # the saddlepoint approximation if it gives a probability (it is rough
        # there, but far better than the failed value). Either way, the value
        # is still flagged by the warning below, since neither has been
        # checked against the other.
        failed <- which(not_conv & !not_rel & !near_one & !small)
        for (k in failed) {
          j  <- idx[k]
          ti <- ptweedie_tail_integrate(q          = q[j],
                                        mu         = mu[j],
                                        phi        = phi[j],
                                        power      = power,
                                        lower.tail = lower.tail,
                                        fallback   = TRUE)
          if (ti$ok && is.finite(ti$logp)) {
            f[j] <- if (log.p) ti$logp else exp(ti$logp)
            next
          }
          f_sp <- ptweedie_saddle(q          = q[j],
                                  mu         = mu[j],
                                  phi        = phi[j],
                                  power      = power,
                                  lower.tail = lower.tail,
                                  log.p      = log.p)
          is_prob <- !is.na(f_sp) && (if (log.p) f_sp < 0 else (f_sp > 0 && f_sp < 1))
          if (is_prob) f[j] <- f_sp
        }

        # A tail probability below 1e-10 that could not be integrated: the
        # inversion's value is accurate only in absolute terms (about 1e-16),
        # so below about 1e-14 it has no correct digits, and can be zero or
        # even negative. There the saddlepoint approximation, where it gives a
        # probability, has at least the right order of magnitude, so it is
        # used instead (still flagged below).
        for (k in which(not_rel & (is.na(tail_lin) | tail_lin < S$p_tail_nodigits))) {
          j    <- idx[k]
          f_sp <- ptweedie_saddle(q          = q[j],
                                  mu         = mu[j],
                                  phi        = phi[j],
                                  power      = power,
                                  lower.tail = lower.tail,
                                  log.p      = log.p)
          is_prob <- !is.na(f_sp) && (if (log.p) f_sp < 0 else (f_sp > 0 && f_sp < 1))
          if (is_prob) f[j] <- f_sp
        }

        n_bad <- sum(not_conv & !not_rel)
        if (n_bad > 0) {
          warning("ptweedie: the numerical integration did not reach the target accuracy for ", n_bad,
                  " value(s); these may be inaccurate.", call. = FALSE)
        }
        n_rel <- sum(not_rel)
        if (n_rel > 0) {
          warning("ptweedie: ", n_rel, " small tail probability(ies) could not be computed to ",
                  "full relative accuracy; these are only approximate (from the saddlepoint ",
                  "approximation, or accurate only in absolute terms).",
                  call. = FALSE)
        }
      }
    } else {
      # CASE 1 < p < 2: the series, at every lambda; the inversion only where
      # the series returns no value.
      # (The choice before 2026 -- by p, after bug reports from Gustavo
      # Lacerda in 2017 -- is in the git history.)
      
      # Since 2026 the SERIES is used first. For 1 < p < 2 it sums positive
      # terms (no cancellation), and it sums exactly the terms that matter
      # (their logs are concave in N, so the peak is found and the sum
      # extended outwards from it), so it has full relative accuracy wherever
      # it can be used. The inversion, by contrast, has only absolute
      # accuracy; it can be silently wrong (errors to about 1e-2) when the
      # point mass at zero is large (small lambda) and p is near 1, and it
      # degrades as lambda grows (relative errors near 5e-8 at lambda = 1e7,
      # and some values of 0.5 in place of the answer for lambda >= 1e7).
      # The number of terms grows only like sqrt(lambda): about 15 sqrt(lambda)
      # at y = mu, at most about 18 sqrt(lambda), and fewer in the tails;
      # about 0.4 microseconds each (0.2 s per value at lambda = 1e9, 4 s at
      # lambda = 1e12). The inversion is no alternative for large lambda: for
      # lambda >= 1e10 it fails (returning about 0.5) for y below the mean. So
      # the series is used for every lambda, with its cap on the number of
      # terms (ptweedie_series(max.terms = )) set from lambda: 25 sqrt(lambda),
      # and at least 1e6.
      if ( any(!special_y_cases)) {
        idx     <- which(!special_y_cases)
        lambda  <- mu[idx]^(2 - power) / (phi[idx] * (2 - power))
        use_ser <- rep(TRUE, length(idx))

        if (any(use_ser)) {
          if (verbose) cat("- With 1 < p < 2: using the series for", sum(use_ser), "value(s)\n")
          j <- idx[use_ser]
          max_terms <- max(1e6, ceiling(25 * sqrt(max(lambda[use_ser]))))
          f[j] <- as.numeric( ptweedie_series(q          = q[j],
                                              mu         = mu[j],
                                              phi        = phi[j],
                                              power      = power,
                                              lower.tail = lower.tail,
                                              log.p      = log.p,
                                              max.terms  = max_terms) )
          # In case the series could not be evaluated, fall back to the inversion
          use_ser[use_ser] <- !is.na(f[j])
        }

        use_inv <- !use_ser
        if (any(use_inv)) {
          if (verbose) cat("- With 1 < p < 2: using the inversion for", sum(use_inv), "value(s)\n")
          j <- idx[use_inv]
          f_TMP <- suppressWarnings(
                     ptweedie_inversion(q          = q[j],
                                        mu         = mu[j],
                                        phi        = phi[j],
                                        power      = power,
                                        lower.tail = lower.tail,
                                        log.p      = log.p,
                                        verbose    = verbose,
                                        details    = TRUE) )
          f[j] <- f_TMP$cdf

          failed <- is.na(f_TMP$cdf) | (f_TMP$exitstatus == 1L)
          n_bad  <- sum(failed)
          if (n_bad > 0) {
            warning("ptweedie: the numerical integration did not reach the target accuracy for ", n_bad,
                    " value(s); these may be inaccurate.", call. = FALSE)
          }
          # (Reached only where the series returned no value.) The
          # inversion's error is ABSOLUTE (and larger than for p > 2: up to
          # about 1e-11 for p near 1), so small tail probabilities have little
          # relative accuracy. For lambda >= 1e10 it also fails (returning
          # about 0.5, flagged by exitstatus) for y below the mean.
          tail_lin <- if (log.p) exp(f_TMP$cdf) else f_TMP$cdf
          n_rel <- sum( !failed & !is.na(tail_lin) & (tail_lin < 1e-5) )
          if (n_rel > 0) {
            warning("ptweedie: ", n_rel, " small tail probability(ies) could not be computed to ",
                    "full relative accuracy; these are accurate only in absolute terms.",
                    call. = FALSE)
          }
        }
      }
    }  
  }
  
  f <- as.vector(f)

  # Sanity fixes: a probability lies in [0, 1], so a log-probability cannot
  # exceed 0. (Rounding in the inversion can give a value fractionally above
  # one, whose log is then fractionally above zero.)
  if (!log.p) {
    f[ !is.na(f) & (f < 0) ] <- 0
    f[ !is.na(f) & (f > 1) ] <- 1
  } else {
    f[ !is.na(f) & (f > 0) ] <- 0
  }

  # The result has the same form (dim, names, ...) as  q, as for pnorm() etc.
  return( match_shape(f, q) )
    
}


################################################################################


#' @export
qtweedie <- function(p, xi = NULL, mu, phi, power = NULL){
  
  ### BEGIN preliminary work
  
  # SORT OUT THE NOTATION (i.e., xi VS power)
  out   <- sort_notation(xi = xi, power = power)
  xi    <- out$xi
  power <- out$power
  
  xi.notation    <- out$xi.notation
  index.par      <- out$index.par
  index.par.long <- out$index.par.long 
  

  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  out <- check_inputs(p, mu, phi, power,
                      type = "quantile")
  mu  <- out$mu
  phi <- out$phi

  # MISSING VALUES: NA in, NA out; the other values are computed as usual
  # (as in ptweedie(); before, an NA gave an error)
  if ( any(is.na(p)) ) {
    ok  <- !is.na(p)
    ans <- rep(NA_real_, length(p))
    if ( any(ok) ) {
      ans[ok] <- qtweedie(p     = p[ok],
                          mu    = mu[ok],
                          phi   = phi[ok],
                          power = power)
    }
    return( match_shape(ans, p) )
  }

  # NO VALUES: a zero-length result (as for qnorm(numeric(0)))
  if ( length(p) == 0 ) return( match_shape(numeric(0), p) )

  f   <- array(0,
               dim = length(p) )
  
  
  # IDENTIFY SPECIAL CASES

  special_y_cases <- rep(FALSE, length(p))
  out <- special_cases(p, mu, phi, power)
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases
  if ( any(special_y_cases) ) {
    special_y_cases <- out$special_y_cases  
    f <- out$f
  }
  
  ### END preliminary work
  
  len <- length(p) 
  
  # Some monkeying around to explicitly account for the cases p=1 and p=0
  ans <- ans2 <- rep( NA, length = len )
  if ( any(p == 1) ) ans2[p == 1] <- Inf
  if ( any(p == 0) ) ans2[p == 0] <- 0
  
  ans     <-  ans[ ( (p > 0) & (p < 1) ) ]
  mu.vec  <-  mu[  ( (p > 0) & (p < 1) ) ]
  phi.vec <-  phi[ ( (p > 0) & (p < 1) ) ]
  p.vec   <-  p[   ( (p > 0) & (p < 1) ) ]
  
  for (i in seq_along(ans) ) {   # (no loop when all p are 0 or 1)
    mu.1  <- mu.vec[i]
    phi.1 <- phi.vec[i]
    p.1   <- p.vec[i]  # This is the  qtweedie()  input p (a probability)
    pwr   <- power     # This is the Tweedie power, xi
    prob  <- p.1 # Rename p to avoid confusion with  pwr: This is the  qtweedie()  input p (a probability)
    
    if ( pwr < 2 ) {
      qp <- stats::qpois(prob, 
                         lambda = mu.1 / phi.1)
      if ( pwr == 1 ) ans[i] <- phi.1 * qp   
    }
    
    qg <- stats::qgamma(prob,  
                        rate  = 1 / (phi.1 * mu.1), 
                        shape = 1 / phi.1 )
    if ( pwr == 2 ) ans[i] <- qg
    
    # Starting values
    # - for 1 < pwr < 2, linearly interpolate between Poisson and gamma
    if ( (pwr > 1) & ( pwr < 2) ) {
      start <- (qg - phi.1 * qp) * pwr + (2 * phi.1 * qp - qg)
    }
    
    # - for pwr > 2, start with gamma
    if ( pwr > 2 ) start <- qg
    
    # Solve!
    if ( ( pwr > 1) & (pwr < 2) ) { # This gives a *lower* bound on the value of the answer (if y>0)
      step <- dtweedie(y     = 0, 
                       mu    = mu.1, 
                       phi   = phi.1, 
                       power = pwr)
      # This is P(Y = 0), the discrete "step"
      
      if ( prob <= step ) {
        ans[i] <- 0
      }
    }
    
    if ( is.na(ans[i]) ) { # All cases except Y = 0 when 1 < pwr < 2
      
      pt2 <- function( q, 
                       mu, 
                       phi, 
                       pwr, 
                       p.given = prob ){ 
        
        ptweedie(q     = q, 
                 mu    = mu, 
                 phi   = phi, 
                 power = pwr ) - p.given
      }
      
      pt <- pt2( q       = start, 
                 mu      = mu.1, 
                 phi     = phi.1, 
                 pwr     = pwr,
                 p.given = prob)
      
      if ( pt == 0 ) {
        # The starting value is the quantile itself.
        # (Before, this was stored in the wrong place, in ans2[i] rather than 
        # ans[i], and uniroot() was then called with  start.2  undefined.)
        ans[i] <- start
        
      } else {
        
        if ( pt > 0 ) { 
          loop <- TRUE
          start.2 <- start
          while ( loop ) {
            # Try harder
            start.2 <- 0.5 * start.2
            if (pt2( q       = start.2, 
                     mu      = mu.1, 
                     phi     = phi.1, 
                     pwr     = pwr, 
                     p.given = prob ) < 0 ) loop <- FALSE
            # RECALL:  We are only in this part of the loop if  pt > 0
          }
        }
        
        if ( pt < 0) {
          loop <- TRUE
          start.2 <- start
          
          while ( loop ) {
            # Try harder
            start.2 <- 1.5 * (start.2 + 2)
            if (pt2( q       = start.2, 
                     mu      = mu.1, 
                     phi     = phi.1, 
                     pwr     = pwr, 
                     p.given = prob ) > 0 ) loop <- FALSE
            # RECALL:  We are only in this part of the loop if  pt < 0
          }
        }
        
        # (The root is bracketed by  start  and  start.2; uniroot() uses the
        # smaller as the lower end of the interval. uniroot() is now called 
        # once only: before, it was also called with its default tolerance,
        # and that result was not used.)
        ans[i] <- stats::uniroot(pt2, 
                                 c(start, start.2), 
                                 mu      = mu.1, 
                                 phi     = phi.1, 
                                 pwr     = pwr, 
                                 p.given = prob, 
                                 tol     = 1e-12 )$root
      }
    }
    
  }
  
  ans2[ is.na(ans2) ] <-  ans

  # The result has the same form (dim, names, ...) as  p, as for qnorm() etc.
  match_shape(ans2, p)
}

################################################################################


#' @export
rtweedie <- function(n, xi = NULL, mu, phi, power = NULL){
  
  ### BEGIN preliminary work

  # SORT OUT THE NOTATION (i.e., xi VS power)
  out <- sort_notation(xi = xi, power = power)
  xi  <- out$xi
  
  power          <- out$power
  xi.notation    <- out$xi.notation
  index.par      <- out$index.par
  index.par.long <- out$index.par.long

  # THE NUMBER OF VALUES: n is used as in R's own r-functions (rnorm(),
  # rgamma(), ...): if  length(n) > 1, the length is the number required;
  # a non-integer  n  is truncated (so n = 2.7 gives two values); and n = 0
  # gives a zero-length result
  if ( length(n) > 1L ) n <- length(n)
  if ( length(n) == 0L ) return( numeric(0) )
  n <- suppressWarnings( as.numeric(n) )
  if ( is.na(n) || !is.finite(n) || (n < 0) ) {
    stop("n must be a non-negative integer.\n")
  }
  n <- floor(n)


  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  out <- check_inputs(y     = n,
                      mu    = mu, 
                      phi   = phi, 
                      power = power,
                      type  = "random")
  mu  <- out$mu
  phi <- out$phi

  if ( n == 0 ) return( numeric(0) )   # (as for rnorm(0))

  # IDENTIFY SPECIAL CASES
  out <- special_cases(n, mu, phi, power)
  f <- out$f
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases
  ### END preliminary work
  
  
  if (power == 1) {
    rtw <- phi * stats::rpois( n, 
                               lambda = mu/phi)
  }
  if (power == 2) {
    alpha <- (2 - power) / (1 - power)
    gam   <- phi * (power - 1) * mu ^ (power - 1)
    rtw   <- stats::rgamma( n, 
                          shape = 1 / phi, 
                          scale = gam )
  }

  if ( power > 2) {
    rtw <- qtweedie( stats::runif(n),
                     mu    = mu,
                     phi   = phi, 
                     power = power)
  }
  
  if ( (power > 1) & (power < 2) ) {
    # Two options:  As above or directly.
    # Directly is faster
    rtw <- numeric(n)
    
    lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
    alpha  <- (2 - power) / (1 - power)
    gam    <- phi * (power - 1) * mu ^ (power - 1)
    
    N <- stats::rpois(n, 
                      lambda = lambda)
    for (i in (1:n) ){
      rtw[i] <- stats::rgamma(1, 
                              shape = -N[i] * alpha, 
                              scale = gam[i])
    }
  }
  # A plain double vector of length n, with no attributes (as for rnorm() etc.)
  as.vector(rtw, mode = "double")
}



