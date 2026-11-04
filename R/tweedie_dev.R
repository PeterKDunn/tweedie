#' @title Unit Deviance for a Tweedie Distribution
#' @name tweedie_dev
#' @description Computes the unit deviance for Tweedie distributions.
#' Since version 3.2.0 it is computed without cancellation, so it has full 
#' relative accuracy even when \code{y} is very close to \code{mu}.
#'
#' @usage tweedie_dev(y, mu, power)
#' @param y vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' 
#' @return A numeric vector containing the unit deviance.
#' 
#' @references
#' Jorgensen, B. (1997).
#' \emph{Theory of Dispersion Models}.
#' Chapman and Hall, London.
#' 
#' @examples
#' # Unit deviance is not symmetric in general:
#' round( tweedie_dev(0:6, mu = 3, power = 1.1), 3)
#' 
#' @export
tweedie_dev <- function(y, mu, power){
  # 
  # Peter K Dunn 
  # 29 Oct 2000 
  # 
  # Rewritten (2026) to avoid cancellation. The unit deviance is small when y
  # is close to mu, but the textbook formula computes it as the difference of
  # terms that are much larger (for 1 < p < 2 and p > 2, of size about
  # 1/{(p - 1)|p - 2|}), so near y = mu it lost most of its digits: a relative
  # error of about 1e-16 / (y/mu - 1)^2, which made exp(-dev/(2 phi)) wrong
  # when phi is small. Instead, with e = log(y/mu) (computed as log1p() of
  # (y - mu)/mu when y is close to mu, where y - mu is exact),
  #
  #   dev = 2 mu^(2 - p) g(e) / {(1 - p)(2 - p)},  g(e) = expm1(b e) - b expm1(e),
  #
  # with a = 1 - p and b = 2 - p. Near p = 1 (a small), g is written instead as
  # exp(e) expm1(a e) - a expm1(e), whose two terms are both of size a. For
  # |e| < 0.1, g still cancels, so its series is summed:
  #
  #   dev = (2 mu^(2 - p) / a) sum_{k >= 2} (b^(k - 1) - 1) e^k / k!,
  #
  # with b^(k - 1) - 1 = expm1((k - 1) log1p(a)) when b > 0. The special cases
  # p = 1 (Poisson) and p = 2 (gamma) are the limits, written the same way.
  # The result has a relative error of about 1e-15 for every y > 0.
  
  p   <- power
  n   <- if ( (length(y) == 0) || (length(mu) == 0) ) 0 else max(length(y), length(mu))
  y   <- rep_len(y, n)
  mu  <- rep_len(mu, n)
  dev <- rep(NA_real_, n)
  
  if ( (p == 0) || (p < 1) ) {
    # The normal distribution (p = 0), and p < 1: the textbook formulas
    if (p == 0) {
      dev <- (y - mu) ^ 2 / 2
    } else{
      dev <- (y ^ (2 - p)) / ((1 - p) * (2 - p)) - 
        (y * (mu ^ (1 - p)))/(1 - p) + 
        (mu ^ (2 - p)) / (2 - p)
    }
    return(dev * 2)
  }
  
  # e = log(y/mu): exactly from y - mu when y is close to mu
  r    <- (y - mu) / mu
  e    <- ifelse( abs(r) < 0.5, 
                  log1p(r), 
                  log(y) - log(mu) )
  pos  <- !is.na(y) & (y > 0)
  near <- pos & is.finite(e) & (abs(e) < 0.1)
  far  <- pos & !near
  
  if (p == 1) {
    # dev/2 = y log(y/mu) - (y - mu) = mu {e exp(e) - expm1(e)}
    dev[far] <- mu[far] * (e[far] * exp(e[far]) - expm1(e[far]))
    if (any(near)) {
      en <- e[near]
      s  <- 0
      tk <- en                                      # e^k / k!
      
      for (k in 2:40) {
        tk <- tk * en / k
        s  <- s + (k - 1) * tk
      }
      dev[near] <- mu[near] * s
    }
    y0      <- !is.na(y) & (y == 0)
    dev[y0] <- mu[y0]
    
    return(dev * 2)
  }
  
  if (p == 2) {
    # dev/2 = y/mu - 1 - log(y/mu) = expm1(e) - e
    dev[far] <- expm1(e[far]) - e[far]
    if (any(near)) {
      en <- e[near]
      s  <- 0
      tk <- en
      
      for (k in 2:40) {
        tk <- tk * en / k
        s  <- s + tk
      }
      dev[near] <- s
    }
    y0      <- !is.na(y) & (y == 0)
    dev[y0] <- Inf
    
    return(dev * 2)
  }
  
  # 1 < p < 2 and p > 2
  a <- 1 - p
  b <- 2 - p
  if ( abs(a) < abs(b) ) {
    g <- exp(e[far]) * expm1(a * e[far]) - a * expm1(e[far])
  } else {
    g <- expm1(b * e[far]) - b * expm1(e[far])
  }
  dev[far] <- mu[far] ^ b * g / (a * b)
  
  if (any(near)) {
    en <- e[near]
    s  <- 0
    tk <- en                                        # e^k / k!
    
    for (k in 2:40) {
      tk <- tk * en / k
      ck <- if (b > 0) expm1( (k - 1) * log1p(a) ) else (b ^ (k - 1) - 1)
      s  <- s + ck * tk
    }
    dev[near] <- mu[near] ^ b * s / a
  }
  
  # Where the formula above overflows (y/mu beyond about exp(700)), or y = 0:
  # the textbook formula, which is then accurate
  redo <- !is.na(y) & ( (y == 0) | (pos & !is.finite(dev)) )
  if (any(redo)) {
    dev[redo] <- (y[redo] ^ (2 - p)) / ((1 - p) * (2 - p)) - 
      (y[redo] * (mu[redo] ^ (1 - p)))/(1 - p) + 
      (mu[redo] ^ (2 - p)) / (2 - p)
  }
  dev * 2
}



#' @rdname tweedie_dev
#' @export
tweedie.dev <- function(y, mu, power){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "tweedie.dev()", 
                            with = "tweedie_dev()")
  tweedie_dev(y, mu, power)
}

