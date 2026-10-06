#' @title Series Evaluation for the Tweedie Distribution Function
#' @name ptweedie_series
#' @description
#' Evaluates the distribution function (\acronym{df}) for Tweedie distributions 
#' with \eqn{1 < p < 2}{1 < p < 2}
#' using an infinite series, for given values of the dependent variable \code{y}, 
#' the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be in the case of evaluation problems.
#'
#' @usage ptweedie_series(q, power, mu, phi, lower.tail = TRUE, log.p = FALSE, verbose = FALSE, details = FALSE)
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
#' @param verbose logical; if \code{TRUE}, displays some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution function and some details.
#' 
#' @return A numeric vector of probabilities (or log-probabilities, if \code{log.p = TRUE}).
#' 
#' @note
#' The 'exact' values for the inverse Gaussian distribution are not really exact, but evaluated using inverse normal distributions,
#' for which very good numerical approximation are available in R.
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
                             lower.tail = lower.tail, log.p = log.p,
                             verbose = verbose, details = details))
    if (details) {
      out <- lapply(names(one[[1]]), function(nm) sapply(one, `[[`, nm))
      names(out) <- names(one[[1]])
      return(out)
    }
    return( vapply(one, as.numeric, numeric(1)) )
  }
  # Here mu and phi are each constant: the summation below expects scalars.
  mu  <- mu[1]
  phi <- phi[1]

    # SET UP
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  drop   <- 39

  # FIND THE LIMITS ON N, the summation index
  # The *lower* limit on N
  lambda_bound <- max(lambda )
  logfmax      <-  -log(lambda_bound)/2
  estlogf      <- logfmax
  N            <- max( lambda_bound )
  
  while ( ( estlogf > (logfmax - drop) ) & ( N > 1 ) ) {
    N <- max(1, N - 2)
    estlogf <- -lambda_bound + N * ( log(lambda_bound) - log(N) + 1 ) - log(N)/2
  }
  lo.N <- max(1, floor(N) )
  
  
  # The *upper* limit on N
  lambda_bound <- min( lambda )
  logfmax      <- -log(lambda_bound) / 2
  estlogf      <- logfmax
  N            <- max( lambda_bound )
  
  while ( estlogf > (logfmax - drop) ) {
    N <- N + 1
    estlogf <- -lambda_bound + N * ( log(lambda_bound) - log(N) + 1 ) - log(N)/2
  }
  hi.N <- max( ceiling(N) )
  if (verbose) cat("Summing over", lo.N, "to", hi.N, "\n")
  
  # ---- Row-wise log-sum-exp helper: never forms the sum in linear space first ----
  rowLogSumExp <- function(M) {
    m <- apply(M, 1, max)
    # guard against a row that is entirely -Inf (e.g. q so small every term underflows)
    m[!is.finite(m)] <- -Inf
    out <- m + log(rowSums(exp(M - m)))
    out[!is.finite(m)] <- -Inf
    out
  }
  
  
  # The Poisson weights alone do not locate the terms that matter: each term
  # is T_N = Pr(N) * G(N), where G is the lower (F_G) or upper (Q_G) incomplete
  # gamma function for the requested tail. Far in the upper tail Q_G grows
  # quickly with N, so the largest T_N lie beyond the largest Poisson weights;
  # likewise, far in the lower tail F_G falls quickly with N, so they lie
  # below them. So start from the Poisson range, then widen it, at either end,
  # until the terms at the ends are negligible (below exp(-drop)) relative to
  # the largest term in each row. The total number of terms is capped.
  max_terms <- 1e6
  capped <- FALSE
  if (hi.N > max_terms) {
    hi.N <- max_terms
    capped <- TRUE
  }
  if (hi.N < lo.N) hi.N <- lo.N

  df_of <- function(N) -2 * alpha * N
  x_vec <- 2 * q / tau
  log_terms <- function(N) {
    # Matrix of log T_N: rows = q, columns = N
    lp <- stats::dpois(N, lambda, log = TRUE)
    lg <- outer(x_vec, df_of(N),
                function(x, df) stats::pchisq(x, df, lower.tail = lower.tail, log.p = TRUE))
    sweep(lg, 2, lp, "+")
  }

  N_lo <- lo.N
  N_hi <- hi.N
  M <- log_terms(N_lo:N_hi)
  # the N = 0 atom, exp(-lambda), belongs to the lower tail only
  log_atom <- if (lower.tail) -lambda else -Inf
  rowmax <- function(M) pmax(apply(M, 1, max), log_atom)
  negligible <- function(col, rm) all( !is.finite(rm) | !(col > rm - drop) )
  # Terms still rising towards the end of the range (in any row) mean the
  # largest terms lie beyond it, even if the end term itself is small (e.g.
  # far in the lower tail with large lambda, where the terms that matter are
  # at small N while the range starts near N = lambda)
  rising <- function(end, nxt) any( is.finite(end) & (end > nxt) )

  # widen upwards
  repeat {
    rm <- rowmax(M)
    k <- ncol(M)
    if ( negligible(M[, k], rm) && !(k > 1 && rising(M[, k], M[, k - 1])) ) break
    if (N_hi >= max_terms) { capped <- TRUE; break }
    step  <- max(10, ceiling((N_hi - N_lo + 1) / 2))
    new_N <- (N_hi + 1):min(max_terms, N_hi + step)
    M     <- cbind(M, log_terms(new_N))
    N_hi  <- max(new_N)
  }
  # widen downwards
  repeat {
    rm <- rowmax(M)
    if ( N_lo <= 1 ) break
    if ( negligible(M[, 1], rm) && !(ncol(M) > 1 && rising(M[, 1], M[, 2])) ) break
    step  <- max(10, ceiling((N_hi - N_lo + 1) / 2))
    new_N <- max(1, N_lo - step):(N_lo - 1)
    M     <- cbind(log_terms(new_N), M)
    N_lo  <- min(new_N)
  }
  if (verbose) cat("Summing over", N_lo, "to", N_hi, "\n")
  if (capped) {
    warning("ptweedie_series: the series was truncated at ", max_terms,
            " terms (lambda = ", signif(lambda, 3), " is too large); ",
            "the result may be inaccurate. Use ptweedie() or ptweedie_inversion() instead.",
            call. = FALSE)
  }
  its <- N_hi - N_lo + 1

  if (lower.tail) {
    # F(y) = exp(-lambda) [the N = 0 atom] + sum_N Pr(N) F_G(y; N)
    M <- cbind(M, rep(log_atom, length(q)))
  }
  # S(y) = sum_N Pr(N) Q_G(y; N): no atom term, since P(Y > y | N = 0) = 0 for y > 0
  logP <- rowLogSumExp(M)
  result <- if (log.p) logP else exp(logP)

  if (details) {
    return( list( cdf = result,
                  iterations = its) )
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
  ptweedie_series(q, power, mu, phi, verbose = FALSE, details = FALSE)
}

