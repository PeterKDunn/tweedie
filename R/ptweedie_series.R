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

  # The summation below shares one range of N across all elements, and is not
  # correct when mu or phi vary between elements (it returned wrong values,
  # with a recycling warning, or stopped with an error). In that case,
  # evaluate element by element.
  n <- max(length(q), length(mu), length(phi))
  if ( (n > 1) && ( (length(unique(mu)) > 1) || (length(unique(phi)) > 1) ) ) {
    q   <- rep_len(q, n)
    mu  <- rep_len(mu, n)
    phi <- rep_len(phi, n)
    one <- lapply(seq_len(n), function(i)
             ptweedie_series(q[i], power = power, mu = mu[i], phi = phi[i],
                             lower.tail = lower.tail, 
                             log.p = log.p,
                             verbose = verbose, 
                             details = details))
    if (details) {
      out <- lapply(names(one[[1]]), 
                    function(nm) sapply(one, `[[`, nm))
      names(out) <- names(one[[1]])
      return(out)
    }
    return( vapply(one, 
                   as.numeric, 
                   numeric(1)) )
  }
  # Here mu and phi are each constant: the summation below expects scalars.
  mu  <- mu[1]
  phi <- phi[1]

  # SET UP
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  drop   <- 39            # terms below exp(-drop) times the largest are negligible
  max_terms <- 1e6        # safety cap on the number of terms summed, per value

  # THE TERMS
  # F(y) = exp(-lambda) + sum_{N >= 1} T_N, with T_N = Pr(N) * F_G(y; N), and
  # S(y) =                sum_{N >= 1} T_N, with T_N = Pr(N) * Q_G(y; N),
  # where F_G and Q_G are the lower and upper incomplete gamma functions,
  # computed through pchisq() (x = 2y/tau on -2 alpha N degrees of freedom).
  #
  # log T_N is CONCAVE in N: log Pr(N) = N log(lambda) - lambda - log(N!) is
  # concave (log-gamma is convex), and so are log F_G and log Q_G as functions
  # of the shape (for integer shapes they are tail probabilities of a Poisson
  # distribution, which is log-concave; checked numerically for general
  # shapes). So the terms have a single peak, which is found directly, and
  # beyond any point past the peak they fall at least geometrically. Hence:
  #   1. find the peak N* (start at N = lambda; double the step until the
  #      terms fall, then bisect on the sign of T_{N+1} - T_N);
  #   2. sum outwards from N* in each direction, stopping at the first term
  #      T_k for which the bound on everything beyond it,
  #        T_k / (1 - r),  r = T_k / T_{k-1} < 1,
  #      is below exp(-drop) times the largest term (and at N = 1 below).
  # Each value of q gets its own range; only the terms that matter are summed.
  # logT(N, x): N and x vectors of the same length (or one of them length 1)
  logT <- function(N, x) {
    stats::dpois(N, lambda, log = TRUE) +
      stats::pchisq(x, -2 * alpha * N, lower.tail = lower.tail, log.p = TRUE)
  }
  log_atom <- if (lower.tail) -lambda else -Inf   # the N = 0 term (lower tail only)

  nq   <- length(q)
  logP <- rep(NA_real_, nq)
  its  <- N_lo <- N_hi <- rep(NA_real_, nq)

  # Values needing no summation
  if (lower.tail) {
    logP[!is.na(q) & (q < 0)]  <- -Inf
    logP[!is.na(q) & (q == 0)] <- -lambda
    logP[!is.na(q) & (q == Inf)] <- 0
  } else {
    logP[!is.na(q) & (q < 0)]  <- 0
    logP[!is.na(q) & (q == 0)] <- log(-expm1(-lambda))
    logP[!is.na(q) & (q == Inf)] <- -Inf
  }
  its[!is.na(logP)] <- 0
  todo <- which(!is.na(q) & (q > 0) & is.finite(q))

  capped <- FALSE
  if (length(todo) > 0) {
    x <- 2 * q[todo] / tau
    m <- length(todo)

    # 1. Peak of each concave sequence: the first N >= 1 with
    #    logT(N + 1) <= logT(N); found for all values at once
    falls <- function(N, xx) !(logT(N + 1, xx) > logT(N, xx))
    N0 <- rep(max(1, round(lambda)), m)
    f0 <- falls(N0, x)
    lo <- hi <- N0
    # peak above N0 (f0 FALSE): gallop up, keeping falls(lo) FALSE
    # peak at/below N0 (f0 TRUE): gallop down, keeping falls(hi) TRUE
    step <- rep(1, m)
    up   <- !f0
    act  <- rep(TRUE, m)
    peak <- rep(NA_real_, m)
    while (any(act)) {
      i <- which(act)
      trial <- ifelse(up[i], hi[i] + step[i], pmax(1, lo[i] - step[i]))
      ft <- falls(trial, x[i])
      # upwards: if it falls there, the bracket is (lo, trial]; else move lo up
      u <- up[i]
      hi[i][u &  ft]  <- trial[u & ft]
      lo[i][u & !ft]  <- trial[u & !ft]
      # downwards: if it does not fall there, the bracket is (trial, hi]; else move hi down
      lo[i][!u & !ft] <- trial[!u & !ft]
      hi[i][!u &  ft] <- trial[!u & ft]
      done_dn1 <- !u & ft & (trial == 1)       # falls already at N = 1: peak is 1
      peak[i[done_dn1]] <- 1
      act[i[(u & ft) | (!u & !ft) | done_dn1]] <- FALSE
      step[i] <- 2 * step[i]
    }
    # bisect where not already settled: falls(lo) FALSE, falls(hi) TRUE
    act <- is.na(peak)
    while (any(act & (hi - lo > 1))) {
      i   <- which(act & (hi - lo > 1))
      mid <- (lo[i] + hi[i]) %/% 2
      fm  <- falls(mid, x[i])
      hi[i][fm]  <- mid[fm]
      lo[i][!fm] <- mid[!fm]
    }
    peak[act] <- hi[act]
    Lmax <- logT(peak, x)

    # 2. Sum outwards from the peak, in chunks, all values at once. Running
    #    log-sum-exp: S = log(sum exp(L - Lmax)) relative to the peak term.
    acc <- rep(1, m)                  # sum of exp(L - Lmax), starting with the peak
    cnt <- rep(1, m)
    nlo <- nhi <- peak
    chunk <- max(16, ceiling(2 * sqrt(lambda)))
    for (dir in c(1, -1)) {
      act  <- is.finite(Lmax)
      if (dir < 0) act <- act & (peak > 1)
      N    <- peak
      prev <- Lmax
      while (any(act)) {
        i  <- which(act)
        off <- seq_len(chunk)
        NN <- outer(N[i], dir * off, "+")                     # rows: values; cols: steps
        valid <- NN >= 1
        NN[!valid] <- 1
        L  <- matrix(logT(as.vector(NN), rep(x[i], chunk)), nrow = length(i))
        L[!valid] <- -Inf
        Lp <- cbind(prev[i], L[, -chunk, drop = FALSE])       # the term before each one
        r  <- exp(L - Lp)
        # bound on the remaining terms from each one on: T_k / (1 - r), on the log scale
        bnd <- ifelse(is.finite(L) & (r < 1), L - log1p(-r), ifelse(is.finite(L), Inf, -Inf))
        stop_here <- (bnd < Lmax[i] - drop) | !valid
        has_stop <- rowSums(stop_here) > 0
        k   <- ifelse(has_stop, max.col(stop_here, ties.method = "first"), chunk)
        use <- (col(L) <= k) & is.finite(L)                  # terms included, per row
        acc[i] <- acc[i] + rowSums(ifelse(use, exp(L - Lmax[i]), 0))
        cnt[i] <- cnt[i] + rowSums(use)
        if (dir > 0) nhi[i] <- N[i] + k else nlo[i] <- pmax(1, N[i] - k)
        N[i]    <- N[i] + dir * chunk
        prev[i] <- L[, chunk]
        over    <- cnt[i] > max_terms
        if (any(over)) capped <- TRUE
        act[i]  <- !has_stop & !over
      }
    }
    lp <- ifelse(is.finite(Lmax), Lmax + log(acc), -Inf)
    # add the atom (lower tail): log(exp(lp) + exp(log_atom))
    if (is.finite(log_atom)) {
      mx <- pmax(lp, log_atom)
      lp <- mx + log(exp(lp - mx) + exp(log_atom - mx))
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
            " terms (lambda = ", signif(lambda, 3), " is too large); ",
            "the result may be inaccurate. Use ptweedie() or ptweedie_inversion() instead.",
            call. = FALSE)
  }
  result <- if (log.p) logP else exp(logP)

  if (details) {
    return( list( cdf        = result,
                  iterations = its,     # number of terms summed (N >= 1), per value
                  N.lo       = N_lo,    # range of N summed, per value
                  N.hi       = N_hi) )
  } else {
    return(result)
  }

}



#' @rdname ptweedie_series
#' @export
ptweedie.series <- function(q, power, mu, phi, verbose = FALSE, details = FALSE){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "ptweedie.series()", 
                            with = "ptweedie_series()")
  ptweedie_series(q, power, mu, phi, verbose = verbose, details = details)
}

