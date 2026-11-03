#' @noRd
ptweedie_saddle <- function(q, mu, phi, power, lower.tail = TRUE, log.p = FALSE) {
  # The Lugannani-Rice saddlepoint approximation to the Tweedie distribution
  # function, as used by ptweedie() for p > 2 when the dispersion is very small
  # (see ptweedie() for where). Vectorised in q, mu and phi.
  #
  #   F(y) ~ Phi(w) + dnorm(w) (1/w - 1/u),
  #   w = sign(y - mu) sqrt(d(y, mu)/phi),   u = t.hat sqrt(phi V(y)),
  #   t.hat = (y^(1-p) - mu^(1-p)) / {phi (1-p)}  (the saddlepoint).
  #
  # Written to avoid cancellation. At small dispersion every y that matters is
  # close to mu, where the textbook formulas subtract nearly equal numbers:
  #
  # - work with e = log(y/mu), after scaling to mean 1 (F(y; mu, phi) =
  #   F(y/mu; 1, phi mu^(p-2))), and use expm1();
  #
  # - the unit deviance d = 2 {expm1(b e) - b expm1(e)}/(a b), with a = 1 - p
  #   and b = 2 - p, still cancels when e is small, so for |e| < 0.1 sum its
  #   series d = 2/(a b) sum_{k >= 2} (b^k - b) e^k / k!;
  #
  # - 1/w - 1/u is 0/0 at y = mu, so for |e| < 1e-5 use its Taylor series
  #   sqrt(phi) {p/6 + p (p - 3) e / 24} (mean 1);
  #
  # - write the small tail as dnorm(w) {Mills(w) + (1/w - 1/u)}, on the log
  #   scale, so it neither cancels nor underflows, with the Mills ratio from
  #   a continued fraction (see mills_ratio()) far in the tails.
  #
  # The relative error (checked against the exact inverse Gaussian, p = 3)
  # falls about 30-fold for each 10-fold increase in xi = mu^(2-p)/{phi (p-2)}:
  # about 1.5e-7 at xi = 1e5, 5e-9 at 1e6, 5e-12 at 1e8, and about 1e-13
  # beyond 1e9. It depends on the dispersion phi mu^(p-2), not on xi alone:
  # near p = 2, xi is large even when the dispersion is not small, and the
  # error is then larger (about 2e-8 at xi = 1e7 when p = 2.001).
  #
  # Returns NA where the approximation breaks down (a tail value <= 0).

  n    <- max(length(q), length(mu), length(phi))
  q    <- rep_len(q, n)
  mu   <- rep_len(mu, n)
  phis <- rep_len(phi * mu^(power - 2), n) # dispersion after scaling to mean 1
  e    <- log(q / mu)
  a    <- 1 - power
  b    <- 2 - power

  d  <- 2 * (expm1(b * e) - b * expm1(e)) / (a * b)
  sm <- !is.na(e) & (abs(e) < 0.1)
  if (any(sm)) {
    es <- e[sm]
    s  <- 0
    tk <- es                              # e^k / k!
    
    for (k in 2:40) {
      tk <- tk * es / k
      s  <- s + (b^k - b) * tk
    }
    d[sm] <- 2 * s / (a * b)
  }

  w    <- sign(e) * sqrt(d / phis)
  u    <- expm1(a * e) / a * exp(power * e / 2) / sqrt(phis)
  corr <- 1 / w - 1 / u
  near <- !is.na(e) & (abs(e) < 1e-5)
  
  if (any(near)) {
    corr[near] <- sqrt(phis[near]) *(power / 6 + power * (power - 3) * e[near] / 24)
  }

  # Lower tail F = Phi(w) + dnorm(w) corr; upper tail S = Phi(-w) - dnorm(w) corr.
  # For the upper tail, change the signs of w and corr and use the same code.
  if (!lower.tail) {
    w    <- -w
    corr <- -corr
  }
  ld    <- stats::dnorm(w, log = TRUE)
  small <- !is.na(w) & (w < 0)            # the requested tail is the small one
  out   <- rep(NA_real_, length(w))

  # Small tail: log{ dnorm(w) (Mills(w) + corr) }, Mills(w) = Phi(w)/dnorm(w)
  arg <- mills_ratio(-w[small]) + corr[small]
  out[small] <- ifelse(arg > 0, 
                       ld[small] + log(pmax(arg, 0)), 
                       NA_real_)

  # Large tail: log{ 1 - dnorm(w) (Mills(-w) - corr) }
  big  <- !is.na(w) & !small
  arg2 <- mills_ratio(w[big]) - corr[big]
  out[big] <- ifelse(arg2 > 0, 
                     log1p(-exp(ld[big] + log(pmax(arg2, 0)))), 
                     NA_real_)

  # y so far out that w overflows: the small tail is 0, the large tail 1
  out[!is.na(w) & (w == -Inf)] <- -Inf
  out[!is.na(w) & (w ==  Inf)] <- 0

  if (log.p) out else exp(out)
}


#' @noRd
mills_ratio <- function(x) {
  # R(x) = Phi(-x) / dnorm(x) for x >= 0, to full relative accuracy.
  # Computed as exp{log Phi(-x) - log dnorm(x)}, the two logs nearly cancel
  # when x is large (both are about -x^2/2), losing about x^2/2 * 1e-16 in
  # absolute terms: 1e-13 at x = 37 (a tail of 1e-300), and everything by
  # x = 1e5 (where the result was even Inf). For x >= 3 use instead the
  # continued fraction R(x) = 1/(x + 1/(x + 2/(x + 3/(x + ...)))), evaluated
  # from the back; 100 terms give full double precision for x >= 3.
  R  <- numeric(length(x))
  lo <- !is.na(x) & (x < 3)
  R[lo] <- exp(stats::pnorm(-x[lo], log.p = TRUE) - stats::dnorm(x[lo], log = TRUE))
  hi <- !is.na(x) & (x >= 3)
  if (any(hi)) {
    xh <- x[hi]
    t  <- 0
    for (k in 100:1) t <- k / (xh + t)
    R[hi] <- 1 / (xh + t)
  }
  R[is.na(x)] <- NA_real_
  R
}
