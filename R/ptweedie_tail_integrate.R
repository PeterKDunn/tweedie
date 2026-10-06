#' @noRd
ptweedie_tail_integrate <- function(q, mu, phi, power, lower.tail = TRUE, IGexact = TRUE) {
  # Tail probability for p > 2 by integrating the density; for scalar q, mu, phi.
  #
  # The Fourier inversion gives F(q) = 1/2 -/+ I/pi, with an error that is
  # controlled in an ABSOLUTE sense (about 1e-15), so tail probabilities much
  # smaller than that have no relative accuracy. The density, however, keeps
  # its relative accuracy far into the tails (it is computed after factoring
  # out exp(-d(y, mu)/(2 phi))). So small tails are computed as
  #   upper: Pr(Y > q)  = f(q) * int_q^Inf f(t)/f(q) dt
  #   lower: Pr(Y <= q) = f(q) * int_0^q   f(t)/f(q) dt
  # Scaling by f(q) makes the integrand O(1) near q, so that integrate()'s
  # relative tolerance is meaningful however small the tail is.
  #
  # Returns a list:
  #   logp: log of the tail probability (on the log scale, so it cannot underflow)
  #   ok:   TRUE if the density and the integration both succeeded
  # (IGexact is passed to dtweedie_inversion; it matters only for p = 3, and
  #  is set to FALSE in the tests to check against the exact inverse Gaussian.)

  dens <- function(t) {
    suppressWarnings(
      dtweedie_inversion(t, mu = mu, phi = phi, power = power,
                         IGexact = IGexact, 
                         details = TRUE) 
      )
  }

  d_q <- dens(q)
  fq  <- d_q$density
  
  if ( !is.finite(fq) || (fq < 0) || (d_q$exitstatus != 0L) ) {
    # The inversion could not compute the density at q. 
    # Where the saddlepoint approximation is accurate (xi = phi q^(p-2) 
    # small; relative error O(xi)) and says the density is far below the
    # smallest double (log density < -800), the density underflows; 
    # and so does the tail probability (e.g. for a lower tail below the 
    # mode, F(q) <= q f(q)). 
    # Report it as 0.
    xi <- phi * q^(power - 2)
    
    if ( is.finite(xi) && (xi < 0.01) ) {
      logf_sad <- -0.5 * log(2 * pi * phi) - (power / 2) * log(q) -
                   tweedie_dev(y = q, 
                               mu = mu, 
                               power = power) / (2 * phi)
      if ( !is.na(logf_sad) && (logf_sad < -800) ) {
        return( list(logp = -Inf, 
                     ok = TRUE) )
      }
    }
    return( list(logp = NA_real_, 
                 ok = FALSE) )
  }
  if ( fq == 0 ) {
    # The density underflows at q, so the tail does too (for practical purposes)
    return( list(logp = -Inf, 
                 ok = TRUE) )
  }

  # UPPER TAIL: Gauss-Laguerre quadrature.
  # Beyond q the density decays roughly like exp(-r (t - q)), where
  #   r = |mu^(1-p) - q^(1-p)| / (phi (p - 1))
  # is the local decay rate of the log density (from the saddlepoint form,
  # log f ~ -d(t, mu)/(2 phi)). With t = q + x/r,
  #   Pr(Y > q) = (1/r) int_0^Inf exp(-x) [exp(x) f(q + x/r)] dx,
  # which is the form Gauss-Laguerre quadrature is built for: n nodes give
  # S ~ (1/r) sum_i w_i exp(x_i) f(q + x_i/r). Rules of 10 and 20 nodes are
  # compared, and the result is accepted only if they agree to 1e-10; this
  # needs 30 density evaluations instead of the several hundred used by the
  # piecewise adaptive integration below, which remains the fallback.
  if (!lower.tail) {
    r_gl <- abs(mu^(1 - power) - q^(1 - power)) / (phi * (power - 1))
    if ( is.finite(r_gl) && (r_gl > 0) ) {
      gl_logS <- function(n) {
        gq <- statmod::gauss.quad(n, kind = "laguerre")
        d  <- dens(q + gq$nodes / r_gl)
        if ( any(d$exitstatus != 0L) || any(!is.finite(d$density)) || any(d$density < 0) ) return(NA_real_)
        lt <- log(gq$weights) + gq$nodes + log(d$density)   # log of each term
        m  <- max(lt)
        if (!is.finite(m)) return(NA_real_)
        m + log(sum(exp(lt - m))) - log(r_gl)
      }
      l10 <- gl_logS(10L)
      l20 <- gl_logS(20L)
      if ( is.finite(l10) && is.finite(l20) && (abs(expm1(l10 - l20)) < 1e-10) ) {
        return( list(logp = l20, ok = TRUE) )
      }
      # otherwise fall through to the adaptive integration
    }
  }

  # The tail mass lies within a few multiples of 1/r of q, where r is the
  # local decay rate of the log density, |d/dy log f|, taken from the
  # saddlepoint form (log f ~ -d(y, mu)/(2 phi)):
  #   r = |mu^(1-p) - q^(1-p)| / (phi (p - 1)).
  # This can be tiny compared with q (e.g. q = 0.01, mu = 0.001, phi = 0.1,
  # p = 2.1: 1/r ~ 6e-5), which integrate() cannot find on [q, Inf) without
  # many thousands of evaluations. 
  #
  # So, integrate over pieces of length
  # (1/r) 2^k moving away from q, until a piece adds a negligible amount.
  # The total number of density evaluations is capped, so that this can never
  # be slow; if the cap is reached, report failure (the caller then keeps the
  # inversion value and warns).
  #
  # Budget: at most 2000 density evaluations and 3 seconds; and stop at the
  # first density value the inversion flags (the result would be rejected).
  max_evals   <- 2000L
  max_secs    <- 3
  t_start     <- proc.time()[["elapsed"]]
  n_evals     <- 0L
  bad_density <- FALSE
  
  g <- function(t) {
    n_evals <<- n_evals + length(t)
    if ( (n_evals > max_evals) ||
         (proc.time()[["elapsed"]] - t_start > max_secs) ) stop("budget exceeded")
    d <- dens(t)
    if (any(d$exitstatus != 0L)) {
      bad_density <<- TRUE
      stop("density not computed accurately")
    }
    d$density / fq
  }

  r <- abs(mu^(1 - power) - q^(1 - power)) / (phi * (power - 1))
  L <- if (is.finite(r) && (r > 0)) {
    1 / r 
  } else {
    q
  }
  if (lower.tail) L <- min(L, q)

  total <- 0
  err   <- 0
  ok    <- TRUE
  k     <- 0L
  repeat {
    a <- (2^k - 1) * L
    b <- (2^(k + 1) - 1) * L
    if (lower.tail) {
      lo <- max(0, q - b)
      hi <- q - a
    } else {
      lo <- q + a
      hi <- q + b
    }
    piece <- tryCatch(stats::integrate(g, 
                                       lower = lo, 
                                       upper = hi,
                                       rel.tol = 1e-10, 
                                       abs.tol = 0, 
                                       subdivisions = 200L),
                      error = function(e) NULL)
    if (is.null(piece)) { 
      ok <- FALSE
      break 
    }
    total <- total + piece$value
    err   <- err + piece$abs.error
    
    if (lower.tail && (lo <= 0)) break                     # reached 0
    if ( (piece$value < 1e-14 * total) && (k >= 2L) ) break # negligible
    k <- k + 1L
    if (k > 200L) { 
      ok <- FALSE
      break 
    }
  }
  out <- list(value = total, 
              abs.error = err)

  if ( !ok || bad_density || 
       !(out$value > 0) ||
       (out$abs.error > 1e-6 * out$value) ) {
    return( list(logp = NA_real_, 
                 ok = FALSE) )
  }

  list(logp = log(fq) + log(out$value), 
       ok = TRUE)
}
