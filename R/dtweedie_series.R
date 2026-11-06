#' @title Series Evaluation for the Tweedie Probability Function
#' @name dtweedie_series
#' @description
#' Evaluates the probability density function (\acronym{pdf}) for Tweedie distributions using an infinite series, 
#' for given values of the dependent variable \code{y}, the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be used in the case of evaluation problems.
#'
#' @usage dtweedie_series(y, power, mu, phi, details = FALSE)
#'
#' @param y vector of quantiles.
#' @param power scalar; the value of \eqn{p}{power} such that the variance is \eqn{\mbox{var}[Y]=\phi\mu^{p}}{var[Y] = phi * mu^power}.
#' @param mu vector of mean \eqn{\mu}{mu}.
#' @param phi vector of dispersion parameters \eqn{\phi}{phi}.
#' @param details logical; if \code{TRUE}, returns a list: the density (\code{density}), 
#'   its logarithm (\code{logf}, computed without underflow), an estimate of its relative error from rounding (\code{relerr}; for \eqn{p > 2}{power > 2}, 
#'   where the series alternates in sign; zero otherwise), and the range of terms summed (\code{lo}, \code{hi}).
#' 
#' @return A numeric vector of probabilities.
#' 
#' @section Accuracy:
#' For \eqn{1 < p < 2}{1 < power < 2}, the terms of the series are all positive, 
#' and the density keeps full relative accuracy.
#' For \eqn{p > 2}{power > 2}, the terms alternate in sign, and rounding error can 
#' cancel many digits. An estimate of the resulting relative error is returned 
#' (as \code{relerr}) when \code{details = TRUE}; no warning is given for a large
#' estimated error. Where the series cannot be summed at all, \code{NaN} is returned 
#' with a warning.
#' \code{\link{dtweedie}()} uses the series for \eqn{p > 2}{power > 2} only where this 
#' estimated relative error is below \eqn{10^{-12}}{1e-12}.
#' 
#' @references
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 

#' @examples 
#' # Plot a Tweedie density
#' y <- seq(0.01, 4, length = 50)
#' fy <- dtweedie_series(y, power = 1.1, mu = 1, phi = 1)
#' plot(y, fy, type = "l", lwd = 2, ylab = "Density")

#' @importFrom stats dgamma dpois 
#' 
#' @keywords distribution
#'
#' @export
dtweedie_series <- function(y, power, mu, phi, details = FALSE){ 
  # Evaluates the Tweedie density using a series expansion
  
  if ( power < 1)      stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y < 0) )    stop("y must be a non-negative vector.")
  if ( any(mu <= 0) )  stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  y0      <- (y == 0)
  yp      <- (y != 0)
  density <- lo <- hi <-array( dim = length(y))
  logf    <- rep(NA_real_, length(y))   # the log of the density
  relerr  <- rep(0, length(y))   # the estimated relative error (p > 2 only; 0 otherwise)
  
  if ( (power == 2) | (power == 1) ) { # Special cases
    if ( power == 2 ){
      density <- dgamma( y, 
                         shape = 1 / phi, 
                         rate  = 1 / (phi * mu ) )
    }
    
    if ( (power == 1) ){
      density <- dpois(x      = y / phi, 
                       lambda = mu / phi )
    }
    # CHANGED: 31 October, thanks to Dina Farkas' email 18 October 2012.  WAS:
    #   if ( (power == 1) & (all(phi==1)) ){
    #    density <- dpois(x=y/phi, lambda=mu/phi )
    #   }
    logf <- log(density)
  } else{
    
    if ( any(y == 0) ) {
      if (power > 2) {
        density[y0] <- 0 * y[y0]
        logf[y0]    <- -Inf
      }
      if ( (power > 1) && (power < 2) ) {
        lambda      <- mu[y0] ^ (2 - power) / ( phi[y0] * (2 - power) )
        density[y0] <- exp( -lambda )
        logf[y0]    <- -lambda
      }
    }
    
    if ( any( y!= 0 ) ) { 
      if (power > 2) {
        out <- dtweedie_series_bigp(power = power,
                                    mu    = mu[yp], 
                                    y     = y[yp],
                                    phi   = phi[yp])
        density[yp] <- out$density 
        logf[yp]    <- out$logf
        relerr[yp]  <- out$relerr
        lo[yp]      <- out$lo
        hi[yp]      <- out$hi
      }
      
      if ( ( power > 1 ) && ( power < 2 ) ) {
        out <- dtweedie_series_smallp(power = power,
                                      mu    = mu[yp], 
                                      y     = y[yp],
                                      phi   = phi[yp])
        density[yp] <- out$density 
        logf[yp]    <- out$logf
        lo[yp]      <- out$lo
        hi[yp]      <- out$hi
      }
    }
  }
  
  # Problematic evaluations:
  n_nan <- sum(is.nan(density))
  if (n_nan > 0) {
    warning("dtweedie_series: the series cannot be summed accurately for ", n_nan,
            " value(s) (catastrophic cancellation, p > 2); NaN returned. ",
            "Use dtweedie() or dtweedie_inversion() for these.", 
            call. = FALSE)
  }
  
  if (details) {
    return(list( density = density,
                 logf    = logf,
                 relerr  = relerr,
                 lo      = lo,
                 hi      = hi) )
  } else {
    return(density)
  }
  
}




#' @rdname dtweedie_series
#' @export
dtweedie.series <- function(y, power, mu, phi){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "dtweedie.series()", 
                            with = "dtweedie_series()")
  dtweedie_series(y, power, mu, phi)
}




#############################################################################


#' @noRd
dtweedie_series_smallp <- function(power, y, mu, phi){ 
  
  #
  # Peter K Dunn
  # 02 Feb 2000
  #
  # Rewritten (2026). For 1 < p < 2, Y is a Poisson sum of gamma variables:
  #
  #   f(y) = sum_{N >= 1} T_N,  T_N = Pr(N; lambda) x g(y; N alpha, tau),
  #
  # with Pr(N; lambda) the Poisson probability (lambda = mu^(2 - p)/{phi (2 - p)})
  # and g the gamma density with shape N alpha (alpha = (2 - p)/(p - 1)) and
  # scale tau = phi (p - 1) mu^(p - 1). This is the same series as before
  # (Dunn and Smyth, 2005), but:
  #
  # - each term is computed as dpois(log = TRUE) + dgamma(log = TRUE), which
  #   R computes without cancellation even for large N (Loader's algorithm);
  #   before, each log term was a difference of numbers of size N log N, and
  #   then the log of the sum had lambda subtracted from it, so the result lost
  #   about log10(lambda) digits (relative errors of about 1e-9 at lambda = 1e4
  #   and 1e-7 at lambda = 1e6);
  #
  # - the terms are summed for each value of y separately, from its own peak
  #   outwards, as in ptweedie_series(). Before, one range of N was used for
  #   all values of y, from the smallest to the largest peak, so a vector y
  #   with widely spread values needed a matrix of (length of y) x (range of
  #   N), which could exhaust the memory.
  #
  # log T_N is concave in N (log Pr(N) is concave, and log g is concave in its
  # shape), so the terms rise to a single peak and then fall: find the peak,
  # then sum outwards from it in each direction until the terms fall below
  # exp(-drop) times the peak term.
  
  if ( power < 1)      stop("power must be between 1 and 2.")
  if ( power > 2)      stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y < 0) )    stop("y must be a non-negative vector.")
  if ( any(mu <= 0) )  stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
  }
  
  drop   <- 39                  # Terms below exp(-drop) times the largest are negligible
  
  n      <- length(y)
  y      <- as.vector(y)            # (plain vectors: 1-d arrays do not recycle)
  lambda <- as.vector( mu ^ (2 - power) / ( phi * (2 - power) ) )
  tau    <- as.vector( phi * (power - 1) * mu ^ ( power - 1 ) )
  alpha  <- (2 - power) / (power - 1)
  
  logf <- rep(NA_real_, n)
  lo   <- hi <- rep(NA_real_, n)
  
  y0 <- !is.na(y) & (y == 0)
  logf[y0] <- -lambda[y0]        # Pr(Y = 0), the point mass
  
  todo <- which( !is.na(y) & (y > 0) & is.finite(y) )
  if (length(todo) > 0) {
    yy  <- y[todo]
    lam <- lambda[todo]
    tt  <- tau[todo]
    m   <- length(todo)
    
    # log T_N, for values i and N (N >= 1)
    logT <- function(N, i) {
      stats::dpois(N, lam[i], log = TRUE) +
        stats::dgamma(yy[i], shape = alpha * N, scale = tt[i], log = TRUE)
    }
    # The derivative of log T_N with respect to N (N continuous); decreasing
    dlogT <- function(N, i) {
      log(lam[i]) - digamma(N + 1) + 
        alpha * ( log(yy[i] / tt[i]) - digamma(alpha * N) )
    }
    
    # 1. The peak: the root of dlogT(N) = 0, found by bisection on log(N),
    #    for all values at once (the peak is N = 1 if the terms fall from N = 1)
    N_lo <- rep(1, m)
    N_hi <- pmax(2, 2 * lam, 2 * yy / (alpha * tt))
    grow <- dlogT(N_hi, seq_len(m)) > 0
    while (any(grow)) {
      N_hi[grow] <- 4 * N_hi[grow]
      grow[grow] <- dlogT(N_hi[grow], which(grow)) > 0
    }
    rising <- dlogT(N_lo, seq_len(m)) > 0
    for (k in 1:60) {
      N_mid <- sqrt(N_lo * N_hi)
      up    <- dlogT(N_mid, seq_len(m)) > 0
      N_lo  <- ifelse(up, 
                      N_mid, 
                      N_lo)
      N_hi  <- ifelse(up, 
                      N_hi, 
                      N_mid)
    }
    N0 <- ifelse(rising, 
                 N_lo, 
                 1)
    
    # The integer peak: the largest of the terms at floor(N0) and floor(N0) + 1
    cand <- cbind( pmax(1, floor(N0)), pmax(1, floor(N0)) + 1 )
    Lc   <- matrix(logT(as.vector(cand), 
                        rep(seq_len(m), 2)), 
                   ncol = 2)
    Lc[is.na(Lc)] <- -Inf
    peak <- cand[cbind(seq_len(m), 
                       max.col(Lc, ties.method = "first"))]
    Lmax <- logT(peak, seq_len(m))
    
    # 2. Sum outwards from the peak, in blocks, relative to the peak term
    acc <- rep(1, m)                         # sum of exp(L - Lmax), starting with the peak
    nlo <- nhi <- peak
    
    for (dir in c(1, -1)) {
      N   <- peak
      act <- which( is.finite(Lmax) & ( (dir > 0) | (peak > 1) ) )
      while (length(act) > 0) {
        B     <- min(1e4, 
                     16 + ceiling(2 * sqrt(max(N[act]))))
        NN    <- outer(N[act], 
                       dir * seq_len(B), 
                       "+")                  # rows: values; cols: steps
        valid      <- NN >= 1
        NN[!valid] <- 1
        L     <- matrix(logT(as.vector(NN), 
                             rep(act, B)), 
                        nrow = length(act)) - Lmax[act]
        keep  <- valid & (L > -drop)         # concavity: once a term is negligible, all beyond are
        nk    <- rowSums(keep)
        acc[act] <- acc[act] + rowSums(exp(L) * keep)
        if (dir > 0) {
          nhi[act] <- N[act] + nk
        } else {
          nlo[act] <- N[act] - nk
        }
        N[act] <- N[act] + dir * B
        act    <- act[ nk == B ]
      }
    }
    logf[todo] <- ifelse(is.finite(Lmax), 
                         Lmax + log(acc), 
                         -Inf)
    lo[todo]   <- nlo
    hi[todo]   <- nhi
  }
  
  list(density = exp(logf), 
       logf    = logf,
       hi      = hi, 
       lo      = lo)
  
}


#############################################################################


#' @noRd
dtweedie_jw_smallp <- function(y, phi, power){ 
  #
  # Peter K Dunn
  # 18 Jun 2002
  #
  
  #
  # Error traps
  #
  
  if ( power < 1)      stop("power must be between 1 and 2.")
  if ( power > 2)      stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be strictly positive.")
  if ( any(y <= 0) )   stop("y must be a strictly positive vector.")
  
  #
  # Set up
  #
  p <- power
  a <- (2 - p) / (1 - p)         # Note that a<0 for 1<p<2
  
  a1 <- 1 - a
  r <- -a * log(y) + a * log(p - 1) - a1 * log(phi) - log(2 - p)        # All terms to power j
  
  drop <- 37             # Accuracy of terms: exp(-37)
  #
  # Find limits of summation using Stirling's approximation
  # to approximate gamma terms, and find max as j.max
  #
  logz  <- max(r)             # To find largest  j  needed
  j.max <- max( y^(2 - p) / ( phi * (2 - p) ) )
  j     <- max( 1, j.max )
  
  c       <- logz + a1 + a * log(-a)
  wmax    <- a1 * j.max
  estlogw <- wmax
  
  # First, the upper limit of j
  
  while(estlogw > (wmax - drop) ){
    j <- j + 2
    estlogw <- j * (c - a1 * log(j))
  }
  
  hi.j <- ceiling(j)
  
  # Now the lower limit of j
  logz  <- min(r) 
  j.max <- min( y ^ ( 2 - power ) / ( phi * (2 - power) ) )
  
  j       <- max( 1, j.max)
  wmax    <- a1 * j.max 
  estlogw <- wmax 
  
  while ( ( estlogw > (wmax - drop) ) && ( j >= 2) ) {
    j       <- max(1, j - 2)
    estlogw <- j*(c - a1 * log(j))
  }
  
  lo.j <- max(1, floor(j))
  # 
  # Now sum the series between established limits.
  # We ensure it works for vector y.
  
  j <- seq( lo.j, hi.j)                 # sequence of j terms needed
  o <- matrix( 1, nrow = length(y))     # matrix of ones
  
  g <- matrix(lgamma( j + 1 ) + lgamma( -a * j ),
              nrow = 1, 
              ncol = hi.j - lo.j + 1) 
  logj <- matrix(log(j),
                 nrow = 1, 
                 ncol = hi.j - lo.j + 1) 
  og     <- o %*% g                             # matrix of gamma terms
  ologj  <- o %*% logj 
  A      <- outer(r, j) - og + ologj            # the series, almost ready to sum
  m      <- apply(A, 1, max)                    # avoid overflow; find maximum values
  we     <- exp( A - m )                        # evaluate terms, less max.
  sum.we <- apply( we, 1, sum)                  # sum terms
  jw     <- sum.we * exp( m )                   # now restore max.
  # Since derivs may be negative, can't use log-scale
  
  list(lo    = lo.j, 
       hi    = hi.j, 
       jw    = jw, 
       j.max = j.max )
  
}



#############################################################################


#' @noRd
dtweedie_kv_bigp <- function(y, phi, power){ 
  # 
  # Peter K Dunn 
  # 18 Jun 2002
  # 
  
  #
  # Error traps
  #
  
  if ( power < 2)      stop("power must be greater than 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) )   stop("y must be a strictly positive vector.")
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  #
  # Set up
  #
  
  p  <- power
  a  <- (2 - p) / (1 - p)
  a1 <- 1 - a
  
  # For p > 2 the series alternates in sign. Its terms peak near
  # k = y^(2-p) / (phi (p-2)), at a size of about exp(k (1-a)), while the sum
  # can be far smaller, so the sum suffers cancellation. Once that peak index
  # is large, the cancellation is beyond double precision (and for very large
  # values the search for the summation limits below never terminates, since
  # k + 2 == k). Values that cannot be summed reliably are returned as NaN.
  # The cancellation is about exp{k_max/(p - 1)} relative to the result, so
  # beyond k_max of about 23 (p - 1) double precision cannot give six correct
  # digits: only values with k_max below 40 (p - 1) + 10 are tried (before, up
  # to 1e6; since one range of k serves all values of y, a single such value
  # made the matrices below enormous).
  n_all     <- length(y)
  kv_all    <- rep(NaN, n_all)
  kmax_each <- y ^ (2 - p) / ( phi * (p - 2) )
  can_try   <- is.finite(kmax_each) & (kmax_each <= 40 * (p - 1) + 10)
  if ( !any(can_try) ) {
    return( list(lo    = NA, 
                 hi    = NA, 
                 kv    = kv_all, 
                 k.max = NA) )
  }
  y   <- y[can_try]
  phi <- phi[can_try]
  
  r    <- -a1 * log(phi) - log(p - 2) - a * log(y) + a * log(p - 1)
  drop <- 37
  
  #
  # Now we find limits of summation, using Stirling's approximation
  # to approximate the gamma terms, and find max as k.max.
  #
  
  logz  <- max(r)
  k.max <- max( y ^ (2 - p) / ( phi * (p - 2) ) )
  k     <- max( 1, k.max )
  
  c       <- logz + a1 + a * log(a)
  vmax    <- k.max * a1
  estlogv <- vmax
  #
  # Now we search either side for when we can
  # ignore terms as they are negligible.
  # Remember we are using the envelope as the tests,
  # not individual terms in the series.
  
  # First:  the upper limit of k
  
  while ( estlogv > (vmax - drop) ) {
    k       <- k + 2
    estlogv <- k * ( c - a1 * log(k) )
  }
  
  hi.k <- ceiling(k)
  #
  # Now the lower limit of k
  #
  logz    <- min(r)
  k.max   <- min( y ^ (2 - p) / ( phi * (p - 2) ) ) 
  k       <- max( 1, k.max )
  c       <- logz + a1 + a * log(a)
  vmax    <- k.max * a1
  estlogv <- vmax
  
  while ( (estlogv > (vmax - drop) ) && ( k >= 2) ) {
    k       <- max(1, k - 2)
    estlogv <- k * ( c - a1 * log(k) )
  }
  
  lo.k <- max(1, floor(k) )
  
  #
  # Now sum the series  between established limits.
  # We ensure it works for vector y.
  #
  k <- seq(lo.k, hi.k) 
  o <- matrix( 1, nrow = length(y)) 
  
  g <- matrix( lgamma( 1 + a * k) - lgamma(1 + k),
               nrow = 1, 
               ncol  =length(k) )
  logk <- matrix( log(k),
                  nrow = 1, 
                  ncol = length(k) )
  
  og    <- o %*% g
  ologk <- o %*% logk
  
  A <- outer(r, k) + og + ologk
  
  C <- matrix( sin( -a * pi * k ) * (-1) ^ k,
               nrow = 1, 
               ncol = length(k) )
  C <- o %*% C
  
  m      <- apply(A, 1, max)
  ve     <- exp(A - m)
  sum.ve <- apply( ve*C, 1, sum )
  # Estimate the rounding error in the sum: each term carries a relative
  # error of about eps * (|A| + 1) (from forming it via exp(A)), and summing
  # n terms adds about eps * n times the sum of the absolute terms. If that
  # error is not small compared to the sum itself, the result is noise
  # (previously returned as 0, a huge value, or Inf): give NaN.
  # (Only the terms that are not negligible count: one range of k serves
  # all values of y, and the largest |A| over that range, or its length,
  # could otherwise swamp the estimate for a value whose terms are few.)
  abs.term        <- abs(ve * C)
  n.eff           <- rowSums(ve > .Machine$double.eps)
  err.est         <- .Machine$double.eps * ( rowSums( (abs(A) + 1) * abs.term ) +
                                               n.eff * rowSums(abs.term) )
  reliable        <- (err.est < 1e-6 * abs(sum.ve))
  kv              <- rep(NaN, length(sum.ve))
  kv[reliable]    <- sum.ve[reliable] * exp( m[reliable] )
  # Since derivs may be negative, can't use log-scale
  kv_all[can_try] <- kv
  
  list(lo    = lo.k, 
       hi    = hi.k, 
       kv    = kv_all, 
       k.max = k.max )
  
}
#############################################################################

#' @noRd
dtweedie_logv_bigp <- function( y, phi, power){ 
  # Peter K Dunn 
  # 02 Feb 2000 
  # 
  
  #
  # Error traps
  #
  
  if ( power < 2)      stop("power must be greater than 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) )   stop("y must be a strictly positive vector.")
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  #
  # Set up
  #
  
  p  <- power
  a  <- (2 - p) / (1 - p)
  a1 <- 1 - a
  
  # For p > 2 the series alternates in sign. Its terms peak near
  # k = y^(2-p) / (phi (p-2)), at a size of about exp(k (1-a)), while the sum
  # can be far smaller, so the sum suffers cancellation. Once that peak index
  # is large, the cancellation is beyond double precision (and for very large
  # values the search for the summation limits below never terminates, since
  # k + 2 == k). Values that cannot be summed reliably are returned as NaN.
  # The cancellation is about exp{k_max/(p - 1)} relative to the result, so
  # beyond k_max of about 23 (p - 1) double precision cannot give six correct
  # digits: only values with k_max below 40 (p - 1) + 10 are tried (before, up
  # to 1e6; since one range of k serves all values of y, a single such value
  # made the matrices below enormous).
  n_all     <- length(y)
  logv_all  <- rep(NaN, n_all)
  kmax_each <- y ^ (2 - p) / ( phi * (p - 2) )
  can_try   <- is.finite(kmax_each) & (kmax_each <= 40 * (p - 1) + 10)
  if ( !any(can_try) ) {
    return( list(lo = NA, 
                 hi = NA, 
                 logv = logv_all, 
                 relerr = rep(NaN, n_all)) )
  }
  y   <- y[can_try]
  phi <- phi[can_try]
  
  r    <- -a1 * log(phi) - log(p - 2) - a * log(y) + a * log(p - 1)
  drop <- 37
  
  #
  # Now we find limits of summation, using Stirling's approximation
  # to approximate the gamma terms, and find max as k.max.
  #
  
  logz  <- max(r)
  k.max <- max( y ^ (2 - p) / ( phi * (p - 2) ) )
  k     <- max( 1, k.max )
  
  c       <- logz + a1 + a * log(a)
  vmax    <- k.max * a1
  estlogv <- vmax
  #
  # Now we search either side for when we can
  # ignore terms as they are negligible.
  # Remember we are using the envelope as the tests,
  # not individual terms in the series.
  
  # First:  the upper limit of k
  
  while ( estlogv > (vmax - drop) ) {
    k       <- k + 2
    estlogv <- k*( c - a1 * log(k) )
  }
  
  hi.k <- ceiling(k)
  #
  # Now the lower limit of k
  #
  logz    <- min(r)
  k.max   <- min( y ^ (2 - p) / ( phi * (p - 2) ) )
  k       <- max( 1, k.max )
  c       <- logz + a1 + a * log(a)
  vmax    <- k.max * a1
  estlogv <- vmax
  
  while ( (estlogv > (vmax - drop) ) && ( k >= 2) ) {
    k       <- max(1, k - 2)
    estlogv <- k*( c - a1 * log(k) )
  }
  
  lo.k <- max(1, floor(k) )
  
  #
  # Now sum the series between established limits.
  # We ensure it works for vector y.
  #
  k <- seq(lo.k, hi.k) 
  o <- matrix( 1, nrow = length(y)) 
  
  g <- matrix( lgamma( 1 + a * k) - lgamma(1 + k),
               nrow = 1, 
               ncol = length(k) )
  
  og <- o %*% g
  A  <- outer(r, k) + og
  C  <- matrix( sin( -a * pi * k ) * (-1)^k,
                nrow = 1, 
                ncol = length(k) )
  
  C      <- o %*% C
  m      <- apply(A, 1, max)
  ve     <- exp(A - m)
  sum.ve <- apply( ve * C, 1, sum )
  
  # Now be careful!  Because of the +/- nature of the sin term,
  # sum.ve can be very small (or even negative) due to subtractive
  # cancellation.
  # Estimate the rounding error in the sum: each term carries a relative
  # error of about eps * (|A| + 1) (from forming it via exp(A)), and summing
  # n terms adds about eps * n times the sum of the absolute terms. If that
  # error is not small compared to the sum itself, the result is noise
  # (previously returned as 0, a huge value, or Inf): give NaN.
  # (Only the terms that are not negligible count: one range of k serves
  # all values of y, and the largest |A| over that range, or its length,
  # could otherwise swamp the estimate for a value whose terms are few.)
  abs.term         <- abs(ve * C)
  n.eff            <- rowSums(ve > .Machine$double.eps)
  err.est          <- .Machine$double.eps * ( rowSums( (abs(A) + 1) * abs.term ) +
                                                n.eff * rowSums(abs.term) )
  reliable         <- (sum.ve > 0) & (err.est < 1e-6 * sum.ve)
  logv             <- rep(NaN, length(sum.ve))
  logv[reliable]   <- log( sum.ve[reliable] ) + m[reliable]
  logv_all[can_try] <- logv
  relerr_all        <- rep(NaN, n_all)      # the estimated relative error of the sum
  relerr_all[can_try][reliable] <- err.est[reliable] / sum.ve[reliable]
  
  list(lo     = lo.k, 
       hi     = hi.k, 
       logv   = logv_all, 
       relerr = relerr_all,
       k.max  = k.max )
  
}


#############################################################################

#' @noRd
dtweedie_logw_smallp <- function(y, phi, power){ 
  #
  # Peter K Dunn
  # 02 Feb 2000
  #
  
  #
  # Error traps
  #
  
  if ( power < 1)      stop("power must be between 1 and 2.")
  if ( power > 2)      stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) )   stop("y must be a strictly positive vector.")
  
  #
  # Set up
  #
  p <- power
  a <- (2 - p) / (1 - p)         # Note that a<0 for 1<p<2
  
  a1 <- 1 - a
  r  <- -a * log(y) + a * log(p - 1) - a1 * log(phi) - log(2 - p)        # All terms to power j
  
  drop <- 37             # Accuracy of terms: exp(-37)
  #
  # Find limits of summation using Stirling's approximation
  # to approximate gamma terms, and find max as j.max
  #
  logz  <- max(r)             # To find largest  j  needed
  j.max <- max( y ^ (2 - p) / ( phi * (2 - p) ) )
  j     <- max( 1, j.max )
  
  cc <- logz + a1 + a * log(-a)    #
  # This is all the terms to j-power, not needing any other j terms.
  # The other terms are introduced when we know the values of j
  
  wmax    <- a1 * j.max
  estlogw <- wmax
  
  # First, the upper limit of j
  while(estlogw > (wmax - drop) ){
    j       <- j + 2
    estlogw <- j * (cc - a1 * log(j))
  }
  
  hi.j <- ceiling(j)
  
  #
  #
  # Now the lower limit of j
  #
  #
  logz  <- min(r) 
  j.max <- min( y ^ (2 - power) / ( phi * (2 - power) ) )
  
  j       <- max( 1, j.max)
  wmax    <- a1 * j.max 
  estlogw <- wmax 
  
  # First, optimize to find the location of the maximum
  while ( ( estlogw > (wmax - drop) ) && ( j >= 2) ) {
    j          <- max(1, j - 2)
    oldestlogw <- estlogw
    estlogw    <- j * (cc - a1 * log(j))
  }
  
  lo.j <- max(1, floor(j))
  # 
  # Now sum the series between established limits.
  # We ensure it works for vector y.
  
  j <- seq( lo.j, hi.j)      # sequence of j terms needed
  
  o <- matrix( 1, nrow = length(y))     # matrix of ones
  
  g <- matrix(lgamma( j + 1 ) + lgamma( -a * j ),
              nrow = 1, 
              ncol = hi.j - lo.j + 1)
  
  og     <- o %*% g              # matrix of gamma terms
  A      <- outer(r, j) - og     # the series, almost ready to sum
  m      <- apply(A, 1, max)     # avoid overflow; find maximum values
  we     <- exp( A - m )         # evaluate terms, less max.
  sum.we <- apply( we, 1, sum)   # sum terms
  logw   <- log( sum.we ) + m    # now restore max.
  
  list(lo    = lo.j, 
       hi    = hi.j, 
       logw  = logw, 
       j.max = j.max )
  
}



#############################################################################

#' @noRd
dtweedie_series_bigp <- function(power, y, mu, phi){ 
  
  # 
  # Peter K Dunn 
  # 02 Feb 2000 
  # 
  
  #
  # Error traps
  #
  
  if ( power < 2)      stop("power must be greater than 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) )   stop("y must be a strictly positive vector.")
  if ( any(mu <= 0) )  stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  result <- dtweedie_logv_bigp(power = power, 
                               y     = y, 
                               phi   = phi)
  logv <- result$logv
  
  theta <- mu ^ (1 - power) / ( 1 - power )
  kappa <- mu ^ (2 - power) / ( 2 - power )
  
  logfnew <- (y * theta - kappa) / phi - log( pi * y) + logv
  f       <- exp( logfnew )
  
  return( list(density = f,
               logf    = logfnew,
               logv    = logv, 
               relerr  = result$relerr,
               lo      = result$lo, 
               hi      = result$hi ) )
  
}
