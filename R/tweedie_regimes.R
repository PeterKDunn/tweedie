################################################################################
### THE REGIMES: which method dtweedie() and ptweedie() use, and where
###
### All the thresholds are in tweedie_regime_settings(), in one place, so they
### can be changed here without touching dtweedie() or ptweedie(). The
### functions dtweedie_regime() and ptweedie_regime() say which regime each
### value falls in; dtweedie() and ptweedie() use them, and they can be called
### directly to see where each value goes:
###
###    tweedie:::dtweedie_regime(y, mu, phi, power)
###    tweedie:::ptweedie_regime(q, mu, phi, power)
###
### The evidence for each threshold is in the comments in dtweedie() and
### ptweedie() (and in the paper, for ptweedie()).
################################################################################


#' @noRd
tweedie_regime_settings <- function() {
  list(
    ### dtweedie()
    # Use the corrected saddlepoint approximation where its estimated relative
    # error is below this (its estimate depends only on p and z = phi y^(p-2))
    d_tol_saddle  = 1e-13,
    # For p > 2, use the series where its estimated relative error (from
    # rounding; the terms alternate in sign) is below this
    d_tol_series  = 1e-12,

    ### ptweedie(), p > 2; xi = mu^(2 - p)/{phi (p - 2)}, phis = phi mu^(p - 2)
    # xi above this: the dispersion is small, and the inversion is not trusted
    # in the tails
    p_xi_trust    = 1e4,
    # xi <= p_xi_trust: tail probabilities (from the inversion) below this
    # are recomputed by integrating the density...
    p_tail_small  = 1e-5,
    # ... and, between p_tail_tiny and p_tail_small, the integration is kept
    # only if it agrees with the inversion to p_agree
    p_tail_tiny   = 1e-10,
    p_agree       = 1e-12,
    # xi > p_xi_trust: tail probabilities below this are integrated
    p_tail_trust  = 1e-2,
    # xi above this: the body is integrated too
    p_xi_body_int = 1e5,
    # xi above this, with phis <= p_phis_saddle: the body uses the saddlepoint
    # approximation
    p_xi_body_sp  = 1e6,
    p_phis_saddle = 1e-5,
    # xi above this (or equal to it, with phis <= p_phis_saddle): the
    # saddlepoint approximation for every value
    p_xi_saddle   = 1e7,
    # A failed tail probability below this (the inversion's value then has
    # no correct digits) is replaced by the saddlepoint approximation
    p_tail_nodigits = 1e-14
  )
}



#' @noRd
dtweedie_regime <- function(y, mu, phi, power) {
  # For each value: the method dtweedie() uses
  #   "exact"         special values of p (1, 2, 3), or y <= 0, or y = Inf
  #   "saddlepoint"   the corrected saddlepoint approximation
  #   "series"        the series (1 < p < 2)
  #   "series/inversion"  p > 2: the series where its estimated error is
  #                   below d_tol_series, otherwise the inversion (which is
  #                   known only once the series has been summed)
  # NA where y is NA.
  S   <- tweedie_regime_settings()
  n   <- max(length(y), length(mu), length(phi))
  y   <- rep_len(y,   n)
  mu  <- rep_len(mu,  n)
  phi <- rep_len(phi, n)

  regime <- rep(NA_character_, n)
  if ( power %in% c(0, 1, 2, 3) ) {
    regime[!is.na(y)] <- "exact"
    return(regime)
  }
  regime[ !is.na(y) & ( (y <= 0) | (y == Inf) ) ] <- "exact"

  pos <- which( !is.na(y) & (y > 0) & (y < Inf) )
  if ( length(pos) > 0 ) {
    sc <- saddle_correction(z = phi[pos] * y[pos] ^ (power - 2),
                            power = power)
    use_sp <- !is.na(sc$relerr) & (sc$relerr < S$d_tol_saddle)
    regime[pos] <- ifelse(use_sp,
                          "saddlepoint",
                          if (power < 2) "series" else "series/inversion")
  }
  regime
}



#' @noRd
ptweedie_regime <- function(q, mu, phi, power) {
  # For each value: the regime ptweedie() uses (as in Table 4 of the paper)
  #   "exact"         special values of p (1, 2, 3), or q <= 0
  #   "series"        1 < p < 2: the series (the inversion where it fails)
  #   For p > 2, by xi = mu^(2 - p)/{phi (p - 2)}:
  #   "inversion"     xi <= p_xi_trust: the inversion; tails below
  #                   p_tail_small by integrating the density
  #   "integration in tails"  p_xi_trust < xi <= p_xi_body_int: tails below
  #                   p_tail_trust by integrating the density; the inversion
  #                   in the body
  #   "integration"   p_xi_body_int < xi <= p_xi_saddle: the density is
  #                   integrated (in the body, the saddlepoint approximation
  #                   when xi > p_xi_body_sp and phis <= p_phis_saddle)
  #   "saddlepoint"   xi > p_xi_saddle (or xi = p_xi_saddle, to within
  #                   rounding, when phis <= p_phis_saddle)
  # Which tail a value is in is decided by the inversion's value, so is
  # known only once that is computed. NA where q is NA.
  S    <- tweedie_regime_settings()
  n    <- max(length(q), length(mu), length(phi))
  q    <- rep_len(q,   n)
  mu   <- rep_len(mu,  n)
  phi  <- rep_len(phi, n)

  regime <- rep(NA_character_, n)
  if ( power %in% c(0, 1, 2, 3) ) {
    regime[!is.na(q)] <- "exact"
    return(regime)
  }
  ok <- !is.na(q)
  if (power < 2) {
    regime[ok] <- "series"
    regime[ok & (q <= 0)] <- "exact"
    return(regime)
  }

  xi   <- mu ^ (2 - power) / (phi * (power - 2))
  phis <- phi * mu ^ (power - 2)
  regime[ok] <- "inversion"
  regime[ok & (xi > S$p_xi_trust)]    <- "integration in tails"
  regime[ok & (xi > S$p_xi_body_int)] <- "integration"
  regime[ok & ptweedie_use_saddle(xi, phis)] <- "saddlepoint"
  regime[ok & (q <= 0)] <- "exact"
  regime
}



#' @noRd
ptweedie_use_saddle <- function(xi, phis) {
  # TRUE where ptweedie() (p > 2) uses the saddlepoint approximation for
  # every value: xi > p_xi_saddle, or xi = p_xi_saddle (to within rounding)
  # with phis <= p_phis_saddle (near p = 2 at xi = p_xi_saddle, the
  # integration is the more accurate)
  S <- tweedie_regime_settings()
  (xi > S$p_xi_saddle * (1 + 1e-8)) |
    ( (xi >= S$p_xi_saddle * (1 - 1e-8)) & (phis <= S$p_phis_saddle) )
}
