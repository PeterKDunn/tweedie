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

    # SET UP
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  drop <- 39

  # FIND THE LIMITS ON N, the summation index
  # The *lower* limit on N
  lambda_bound <- max(lambda )
  logfmax <-  -log(lambda_bound)/2
  estlogf <- logfmax
  N <- max( lambda_bound )
  
  while ( ( estlogf > (logfmax - drop) ) & ( N > 1 ) ) {
    N <- max(1, N - 2)
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  lo.N <- max(1, floor(N) )
  
  
  # The *upper* limit on N
  lambda_bound <- min( lambda )
  logfmax <-  -log(lambda_bound) / 2
  estlogf <- logfmax
  N <- max( lambda_bound )
  
  while ( estlogf > (logfmax - drop) ) {
    N <- N + 1
    estlogf <- -lambda_bound + N * ( log(lambda_bound) - log(N) + 1 ) - log(N)/2
  }
  hi.N <- max( ceiling(N) )
  if (verbose) cat("Summing over", lo.N, "to", hi.N, "\n")
  
  # Add a safety check:
  hi.N <- min(hi.N, 1e6)
  if (hi.N < lo.N) hi.N <- lo.N
  
  N_vec <- lo.N : hi.N
  

  # ---- Row-wise log-sum-exp helper: never forms the sum in linear space first ----
  rowLogSumExp <- function(M) {
    m <- apply(M, 1, max)
    # guard against a row that is entirely -Inf (e.g. q so small every term underflows)
    m[!is.finite(m)] <- -Inf
    out <- m + log(rowSums(exp(M - m)))
    out[!is.finite(m)] <- -Inf
    out
  }
  
  
  # log of the Poisson weights for N = lo.N, ..., hi.N
  log_pois <- dpois(N_vec, 
                    lambda, 
                    log = TRUE)
  
  # log of the incomplete-gamma (via the chi-square relationship already used here),
  # computed directly on the log scale -- never via log(pchisq(...))
  df_vec <- -2 * alpha * N_vec
  x_vec  <- 2 * q / tau
  
  
  # The incomplete Gamma values
  # We want a matrix where rows = q and columns = N
  log_incgamma_lower <- outer(x_vec, 
                              df_vec, 
                              function(x, df) {stats::pchisq(x, df, log.p = TRUE)} )
  log_incgamma_upper <- outer(x_vec, 
                              df_vec, 
                              function(x, df) {stats::pchisq(x, df, lower.tail = FALSE, log.p = TRUE)})
  
  
  # add the (column-wise) log Poisson weights to each column
  log_terms_lower <- sweep(log_incgamma_lower, 
                           2, 
                           log_pois, "+")
  log_terms_upper <- sweep(log_incgamma_upper, 
                           2, 
                           log_pois, "+")
  
  if (!lower.tail) {
    # Upper tail: S(y) = sum_{N=lo.N}^{hi.N} Pr(N) * Q_G(y; N) -- no atom term,
    # since P(Y>y | N=0) = 0 for y>0.
    logS <- rowLogSumExp(log_terms_upper)
    its <- hi.N - lo.N + 1
    result <- if (log.p) logS else exp(logS)
  } else {
    # Lower tail: F(y) = exp(-lambda) [the N=0 atom] + sum_{N=lo.N}^{hi.N} Pr(N) * F_G(y; N)
    log_atom <- -lambda
    # augment with the atom as one extra "term" per row, then logsumexp across all of them
    augmented <- cbind(log_terms_lower, rep(log_atom, length(q)))
    logF <- rowLogSumExp(augmented)
    its <- hi.N - lo.N + 1
    result <- if (log.p) logF else exp(logF)
  }
  
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

