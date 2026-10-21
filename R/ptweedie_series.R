#' @title Series Evaluation for the Tweedie Distribution Function
#' @name ptweedie_series
#' @description
#' Evaluates the distribution function (\acronym{df}) for Tweedie distributions 
#' with \eqn{1 < p < 2}{1 < p < 2}
#' using an infinite series, for given values of the dependent variable \code{y}, 
#' the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be in the case of evaluation problems.
#'
#' @usage ptweedie_series(q, power, mu, phi, lower.tail = TRUE, log.p = FALSE, 
#'                        verbose = FALSE, details = FALSE)
#' 
#' @param q vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param lower.tail logical; if \code{TRUE} (the default), computes \eqn{F(y)=\Pr(Y\le y)};
#'   if \code{FALSE}, computes the upper tail \eqn{S(y)=\Pr(Y>y)} directly, rather than as
#'   \code{1 - ptweedie_series(...)}.
#' @param log.p logical; if \code{TRUE}, the logarithm of the probability is returned,
#'   computed directly via log-sum-exp accumulation over the series terms rather than as
#'   \code{log(ptweedie_series(...))}. This retains accuracy even when the probability
#'   itself would underflow to zero.
#' @param verbose logical; if \code{TRUE}, displays some internal computation details. 
#'   The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution function 
#'   and some details.
#' 
#' @return A numeric vector of probabilities (or log-probabilities, if \code{log.p = TRUE}).
#' 
#' @note
#' The 'exact' values for the inverse Gaussian distribution are not really exact, but evaluated using 
#' inverse normal distributions, for which very good numerical approximation are available in R.
#' 
#' @references
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' @examples
#' # Plot a Tweedie distribution function
#' y <- seq(0.01, 4, length = 50)
#' Fy <- ptweedie_series(y, power = 1.1, mu = 1, phi = 1)
#' plot(y, Fy, type = "l", lwd = 2, ylab = "Distribution function")
#' 
#' @importFrom stats dpois pchisq
#'
#' @keywords distribution
#'
#' @export
ptweedie_series <- function(q, power, mu, phi, lower.tail = TRUE, log.p = FALSE, 
                            verbose = FALSE, details = FALSE) {
  ### NOTE: No notation checks
  ### NOTE: Only for 1 < p < 2
  if ( (length(power) != 1) || !is.finite(power) || (power <= 1) || (power >= 2) ) {
    stop("ptweedie_series() requires 1 < power < 2; ",
         "use ptweedie() or ptweedie_inversion() for other values of power.",
         call. = FALSE)
  }
  
  # In the CDF, we do things differently that for the PDF.
  #
  # # THE SERIES
  #   F(y) = exp(-lambda) + sum_{N >= 1} T_N, with T_N = Pr(N) * F_G(y; N), and
  #   S(y) =                sum_{N >= 1} T_N, with T_N = Pr(N) * Q_G(y; N),
  # where Pr(N) is Poisson(lambda) and F_G, Q_G are the lower and upper
  # regularised incomplete gamma functions with shape -alpha*N, computed
  # through pchisq() (x = 2y/tau on -2 alpha N degrees of freedom).
  #
  # We treat N as continuous.
  # log T_N is strictly concave in N; log Pr(N) is concave; and F_G and Q_G 
  # are log-concave in the shape (the density of the log of a gamma variable 
  # is log-concave, so truncating it, from either side, cannot increase its 
  # variance). 
  #
  # So, for each value of q:
  #   1. find the peak N* by solving d log T_N / dN = 0 (a single root);
  #   2. sum outwards from the integer peak in each direction (as for the
  #      pdf series), stopping once a term falls below exp(-drop) times the
  #      peak term; by concavity, every term beyond that point is smaller still.
  #
  # Step 1 uses uniroot() for each value of q in turn.
  # Step 2 is done for all values of q at once.
  
  drop      <- 39         # Terms below exp(-drop) times the largest are negligible
  max_terms <- 1e6        # Safety cap on the number of terms summed, per value
  
  lens <- c(length(q), length(mu), length(phi))
  n    <- if (min(lens) == 0) 0 else max(lens)
  q    <- rep_len(q,   n)
  mu   <- rep_len(mu,  n)
  phi  <- rep_len(phi, n)
  
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  
  logP <- rep(NA_real_, n)
  its  <- N_lo <- N_hi <- rep(NA_real_, n)
  
  # Values needing no summation
  ok <- !is.na(q) & !is.na(lambda)
  
  if (lower.tail) {
    logP[ok & (q < 0)]    <- -Inf
    logP[ok & (q == 0)]   <- -lambda[ok & (q == 0)]
    logP[ok & (q == Inf)] <- 0
  } else {
    logP[ok & (q < 0)]    <- 0
    logP[ok & (q == 0)]   <- log(-expm1(-lambda[ok & (q == 0)]))
    logP[ok & (q == Inf)] <- -Inf
  }
  its[!is.na(logP)] <- 0
  todo <- which(ok & (q > 0) & is.finite(q))
  
  capped <- FALSE
  if (length(todo) > 0) {
    
    m   <- length(todo)
    x   <- 2 * q[todo] / tau[todo]
    lam <- lambda[todo]
    
    # log T_N, and d log T_N / dN (N continuous), for values i: the Poisson
    # part exactly; the incomplete-gamma part by a central difference
    # (pchisq(log.p = TRUE) is accurate far into both tails)
    logT <- function(N, i) {
      stats::dpois(N, 
                   lam[i], 
                   log = TRUE) +
        stats::pchisq(x[i], 
                      -2 * alpha * N, 
                      lower.tail = lower.tail, 
                      log.p = TRUE)
    }
    dlogT <- function(N, i) {
      h <- 1e-4 * N
      log(lam[i]) - digamma(N + 1) +
        ( stats::pchisq(x[i], 
                        -2 * alpha * (N + h),
                        lower.tail = lower.tail, 
                        log.p = TRUE) -
            stats::pchisq(x[i], 
                          -2 * alpha * (N - h), 
                          lower.tail = lower.tail, 
                          log.p = TRUE) ) / (2 * h)
    }
    
    # 1. The peak N*: the root of d log T_N / dN, which is decreasing in N
    Nstar <- rep(1, m)                    # the peak is N = 1 if already falling there
    for (k in seq_len(m)) {
      f   <- function(N) { dlogT(N, k) }
      lo  <- 1
      flo <- f(lo)
      if (flo <= 0) next
      # Bracket: for the lower tail the derivative is negative at N = lambda;
      # for the upper tail the peak can lie well above lambda, so double
      hi  <- max(2, lam[k])
      fhi <- f(hi)
      while (fhi > 0) {
        lo  <- hi
        flo <- fhi
        hi  <- 2 * hi
        fhi <- f(hi)
      }
      Nstar[k] <- stats::uniroot(f, 
                                 c(lo, hi), 
                                 f.lower = flo, 
                                 f.upper = fhi,
                                 tol = 0.25)$root   # only the integer peak is needed
    }
    
    # The integer peak is next to N*: take the largest term nearby
    cand <- outer(floor(Nstar), -1:2, "+")
    cand[cand < 1] <- 1
    Lc   <- matrix(logT(as.vector(cand), 
                        rep(seq_len(m), 4)), 
                   ncol = 4)
    Lc[is.na(Lc)] <- -Inf
    peak <- cand[cbind(seq_len(m), 
                       max.col(Lc, ties.method = "first"))]
    Lmax <- logT(peak, seq_len(m))
    
    # 2. Sum outwards from the peak, in blocks, relative to the term.
    acc   <- rep(1, m)                  # sum of exp(L - Lmax), starting with the peak
    cnt   <- rep(1, m)
    nlo   <- nhi <- peak
    
    for (dir in c(1, -1)) {
      N    <- peak
      act  <- which(is.finite(Lmax) & ( (dir > 0) | (peak > 1) ))
      
      while (length(act) > 0) {
        B     <- min(1e4, 16 + ceiling(2 * sqrt(max(N[act]))))
        NN    <- outer(N[act], 
                       dir * seq_len(B),
                       "+")                     # rows: values; cols: steps
        valid      <- NN >= 1
        NN[!valid] <- 1
        L     <- matrix(logT(as.vector(NN), rep(act, B)), 
                        nrow = length(act)) -	Lmax[act]
        keep  <- valid & (L > -drop)    # concavity: once a term is negligible, all beyond are
        nk    <- rowSums(keep)
        acc[act] <- acc[act] + rowSums(exp(L) * keep)
        cnt[act] <- cnt[act] + nk
        
        if (dir > 0) {
          nhi[act] <- N[act] + nk
        } else {
          nlo[act] <- N[act] - nk
        }
        over    <- cnt[act] > max_terms
        if (any(over)) capped <- TRUE
        N[act] <- N[act] + dir * B
        act    <- act[ (nk == B) & !over ]
      }
    }
    lp <- ifelse(is.finite(Lmax), 
                 Lmax + log(acc), 
                 -Inf)
    # Add the mass at zero (lower tail): log(exp(lp) + exp(log_atom))
    if (lower.tail) {
      mx <- pmax(lp, -lam)
      lp <- mx + log(exp(lp - mx) + exp(-lam - mx))
    }
    # A probability cannot exceed one; rounding in many terms can give log P
    # just above 0 (by about 1e-12) when P is very close to one
    logP[todo] <- pmin(lp, 0)
    its[todo]  <- cnt
    N_lo[todo] <- nlo
    N_hi[todo] <- nhi
  }
  
  if (verbose) cat("Summed", sum(its, na.rm = TRUE), "terms in all\n")
  if (capped) {
    warning("ptweedie_series: the series was truncated at ", max_terms,
            " terms (lambda too large); the result may be inaccurate. ",
            "Use ptweedie() or ptweedie_inversion() instead.",
            call. = FALSE)
  }
  result <- if (log.p) logP else exp(logP)
  
  if (details) {
    return( list( cdf        = result,
                  iterations = its,     # number of terms summed (N >= 1), per value
                  N.lo       = N_lo,    # range of N summed, per value
                  N.hi       = N_hi) )
  }
  
  result
}



#' @rdname ptweedie_series
#' @export
ptweedie.series <- function(q, power, mu, phi, verbose = FALSE, details = FALSE){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "ptweedie.series()", 
                            with = "ptweedie_series()")
  ptweedie_series(q, power, mu, phi, verbose = verbose, details = details)
}