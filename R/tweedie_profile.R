#' @title Profile Likelihood Estimate of Tweedie Variance Index Parameter
#' @name tweedie_profile
#' @description This function profiles the (log-)likelihood over a vector of
#'   Tweedie power-index parameter (denoted \eqn{p}{power} or \eqn{\xi}{xi}) to find the maximum
#'   likelihood estimate (MLE) of the index parameter \eqn{p} (or equivalently \eqn{\xi}{xi}).
#'
#' @usage tweedie_profile(formula, p.vec = NULL, xi.vec = NULL, link.power = 0, 
#'   data, weights = 1, offset = 0, fit.glm = FALSE, do.smooth = TRUE, 
#'   do.plot = FALSE, do.ci = do.smooth, eps = 1/6, 
#'   control = list( epsilon = 1e-09, maxit = stats::glm.control()$maxit, 
#'   trace = glm.control()$trace ),
#'   do.points = do.plot, method = "inversion", conf.level = 0.95, 
#'   phi.method = ifelse(method == "saddlepoint", "saddlepoint", "mle"), 
#'   verbose = FALSE, add0 = FALSE,
#'   plot_args = list(), point_args = list(), line_args = list())
#'   
#' @details
#' For each value in \code{p.vec}, the function computes an estimate of \eqn{\phi}{phi}
#' and then computes the value of the log-likelihood for these parameters.
#' The plot of the log-likelihood against \code{p.vec} allows the maximum
#' likelihood value of \eqn{p}{p} to be found. 
#' Once \eqn{p}{p} is found, the distribution within the class of Tweedie distributions is identified.
#'
#' @note
#' The estimates of \eqn{p}{p} and \eqn{\phi}{phi} are printed invisibly. 
#' If the response variable has any exact zeros, the values in \code{p.vec} must all be
#' between one and two (or perhaps zero).
#'
#' The function can be temperamental (for theoretical reasons involved in numerically computing
#' the density; see Dunn and Smyth (2005)) and may be very slow or fail. 
#' One solution is to change the method. 
#' The default is \code{method = "inversion"}; then try \code{"series"}, \code{"interpolation"}, 
#' and \code{"saddlepoint"} in that order.
#' Note that \code{method = "saddlepoint"} is an approximate method only.
#'
#' It is recommended that for the first use with a data set, use \code{p.vec}
#' with only a small number of values and set \code{do.smooth = FALSE},
#' \code{do.ci = FALSE}. If this is successful, a larger vector \code{p.vec}
#' and smoothing can be used.
#'
#' @param formula a formula expression as for other regression models and
#'   generalized linear models, of the form \code{response ~ predictors}.
#' @param p.vec a vector of \eqn{p}{power} values for consideration. 
#'   The values must all be larger than one. 
#'   If the response has zeros, values must be \eqn{1 < p < 2}.
#'   If \code{NULL} (default), \code{p.vec} is set automatically.
#' @param xi.vec a synonym for \code{p.vec}, as some authors use the \eqn{\xi}{xi} notation.
#' @param link.power the power link function to use in the \code{tweedie} \acronym{glm} family. 
#'  These link functions \eqn{g(\cdot)}{g()} are of the form \eqn{g(\eta)=\eta^{link.power}}{g(\eta) = \eta^{\text{link.power}}}, where \code{link.power = 0} (default) refers to the logarithm link function.
#' @param data an optional data frame, list or environment containing the variables.
#' @param weights an optional vector of prior weights to be used in the fitting process.
#'  As for \code{\link{glm}}, observation \eqn{i} has dispersion \eqn{\phi/w_i}{phi/w[i]};
#'   observations with zero weight are dropped.
#' @param offset an \emph{a priori} known component included in the linear predictor.
#'   See \code{\link{model.offset}}.
#' @param fit.glm logical; if \code{TRUE}, the Tweedie \acronym{glm} is fitted using the
#'   value of \eqn{p} found by the profiling function. The default is \code{FALSE}.
#' @param do.smooth logical; if \code{TRUE} (default), a spline is fitted to the
#'   data to smooth the profile likelihood plot. \bold{Note} that \code{p.vec}
#'   must contain \emph{at least five points} for smoothing.
#' @param do.plot logical; if \code{TRUE}, a plot of the profile likelihood is produced.
#'   The default is \code{FALSE}.
#' @param do.ci logical; if \code{TRUE}, the nominal \code{100*conf.level} is computed.
#'   Defaults to the value of \code{do.smooth}. 
#'   Confidence intervals are only computed if \code{do.smooth = TRUE}.
#' @param eps the offset in computing the variance function. Default is \code{1/6} (as recommended by Nelder and Pregibon, 1987).
#'   \code{eps} is ignored unless \code{method = "saddlepoint"}.
#' @param control a list of parameters for controlling the fitting process;
#' @param do.points logical; if \code{TRUE}, the points used to compute the likelihood as given by \code{p.vec} (or equivalently, \code{xi.vec}) are explicitly shown by points. 
#' The defaults is the value of \code{do.plot}.
#' @param method the method of evaluation; one of \code{saddlepoint}, \code{interpolation}, \code{series} or \code{inversion} (the default).
#'   When \eqn{\phi}{phi} is found by maximum likelihood (\code{phi.method = "mle"}),
#'   the log-likelihood found along the way (using \code{\link{dtweedie}}) is used,
#'   so \code{method} only matters for \code{method = "saddlepoint"}.
#' @param conf.level the level of confidence for the confidence intervals; the default is \code{0.95} (for \eqn{95\%}{95\%} confidence intervals).
#' @param phi.method the method used to estimate \eqn{\phi}{phi}; one of \code{saddlepoint}, \code{mle} (the default).
#' @param verbose logical; if \code{TRUE}, some details of the calculations are shown. The default is \code{FALSE}.
#' @param add0 logical; if \code{TRUE}, adds \eqn{P(Y = 0)}{P(Y = 0)} to the plot. The default is \code{FALSE}.
#' @param plot_args A named list of arguments controlling the main plot.
#'   These are passed to \code{\link[graphics]{plot}}.
#' @param point_args A named list of graphical parameters for plotted points.
#'   These are passed to \code{\link[graphics]{points}}.
#' @param line_args A named list of graphical parameters for plotted lines.
#'   These are passed to \code{\link[graphics]{lines}}.
#' 
#' @references
#' Dunn, P. K. and Smyth, G. K. (2018).
#' Generalized linear models with examples in R.
#' Springer.
#' \doi{10.1007/978-1-4419-0118-7}
#'
#' @examples 
#' data(rock)
#' out <- tweedie_profile(perm~1, data=rock, do.plot=FALSE, 
#'                        xi.vec=seq(1.5, 2.75, length=11),
#'                        line_args = list(lwd = 2))
#'
#' # The estimate for the variance power index (p, or xi) is:
#' out$xi.max
#' 
#' @importFrom graphics lines rug points abline title plot
#' @importFrom stats fitted optimize glm.fit splinefun glm.control uniroot qchisq dgamma dnorm 
#' 
#' @keywords  models
#' 
#' @export
tweedie_profile <- function(formula, 
                            p.vec = NULL, 
                            xi.vec = NULL, 
                            link.power = 0, 
                            data, 
                            weights = 1, 
                            offset = 0, 
                            fit.glm = FALSE, 
                            do.smooth = TRUE, 
                            do.plot = FALSE, 
                            do.ci = do.smooth, 
                            eps = 1/6,
                            control = list( epsilon = 1e-09, 
                                            maxit = stats::glm.control()$maxit, 
                                            trace = glm.control()$trace ),
                            do.points = do.plot, 
                            method = "inversion",
                            conf.level = 0.95, 
                            phi.method = ifelse(method == "saddlepoint", "saddlepoint", "mle"), 
                            verbose = FALSE, 
                            add0 = FALSE,
                            plot_args  = list(),
                            point_args = list(),
                            line_args  = list()) {
  # verbose gives feedback on screen:
  #    0 : minimal (FALSE)
  #    1 : small amount (TRUE)
  #    2 : lots
  # Determine the value of  p  to use to fit a Tweedie distribution
  # A profile likelihood is used (can can be plotted with do.plot=TRUE)
  
  # The plot can be (somewhat) recovered using (where  out  is the
  # returned object, with smoothing):
  #     plot( out$x, out$y)
  # The actual points can then be added:
  #     points( out$p, out$L)
  # (Note that out$L and out$y are the same if smoothing is  not requested.  
  # Likewise for out$p and out$x.)
  
  # Peter Dunn
  # 07 Dec 2004
  
  
  # BEGIN: Check for some given parameters supplied via args
  plot_defaults <- list(
    col = "black"
  )
  point_defaults <- list(pch = 19,
                         col = "black",
                         cex = 1)
  line_defaults  <- list(col = "black", 
                         lwd = 1)
  
  plot_args  <- utils::modifyList(plot_defaults, 
                                  plot_args)
  point_args <- utils::modifyList(point_defaults, 
                                  point_args)
  line_args  <- utils::modifyList(line_defaults, 
                                  line_args)
  # END: Check for some given parameters supplied via args
  
  # Check the method names now: a typo is then an error, rather than
  # (as before) silently using  method = "series".  Partial matching still works.
  # (phi.method's default depends on  method, so  method  is checked first.)
  method     <- match.arg(method, 
                          c("inversion", "series", "interpolation", "saddlepoint"))
  phi.method <- match.arg(phi.method, 
                          c("mle", "saddlepoint"))
  
  
  if ( is.logical( verbose ) ) verbose <- as.numeric(verbose)
  
  if (verbose >= 1 ) {
    cat("---\n This function may take some time to complete.\n")
    cat("If it fails, try using  method=\"series\"\n")
    cat(" rather than the default  method=\"inversion\"\n")
    cat(" Another possible reason for failure is the range of p;\n")
    cat(" try a different input for  p.vec\n---\n")
  }
  
  mf <- match.call()
  m <- match(c("formula", "data", "weights","offset"), names(mf), 0L)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- stats::model.response(mf, "numeric")
  if ( length(dim(Y)) == 1L ) {
    # A 1-d array response (e.g., from tapply()) is made a vector, as glm()
    # does; otherwise glm.fit() fails with "non-conformable arrays"
    nm <- rownames(Y)
    dim(Y) <- NULL
    if ( !is.null(nm) ) names(Y) <- nm
  }
  X <- if (!stats::is.empty.model(mt)) {
    stats::model.matrix(mt, mf, contrasts = NULL)
  } else { 
    matrix(, NROW(Y), 0)
  }
  weights <- as.vector(stats::model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  if (!is.null(offset)) {
    if (length(offset) == 1)
      offset <- rep(offset, NROW(Y))
    else if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  
  # Prior weights: as in glm(), observation i has dispersion  phi / weights[i].
  # Observations with zero weight carry no information, so are dropped
  # (otherwise their dispersion would be  phi / 0 = Inf).
  if ( is.null(weights) ) weights <- rep(1, NROW(Y))
  if ( any(weights == 0) ) {
    keep.obs <- (weights > 0)
    Y        <- Y[keep.obs]
    X        <- X[keep.obs, , drop = FALSE]
    weights  <- weights[keep.obs]
    if ( !is.null(offset) ) offset <- offset[keep.obs]
  }
  
  ### NOW some notation stuff
  if ( !is.null(p.vec) & !is.null(xi.vec) ) {
    warning("Both  p.vec  and  xi.vec  were supplied; only  xi.vec  is used.")
  }
  xi.notation <- TRUE
  if ( is.null(xi.vec) & !is.null(p.vec) ){ # If p.vec given, but not xi.vec
    xi.vec      <- p.vec
    xi.notation <- FALSE
  }

  if ( is.null(xi.vec)) { # Neither xi.vec or p.vec are given
    if ( any(Y == 0) ){
      xi.vec <- seq(1.2, 1.8, 
                    by = 0.1)
    } else {
      xi.vec <- seq(1.5, 5, 
                    by = 0.5)
    }
  }
  
  ### Hereafter, only  xi.vec  is used (whichever notation was supplied)
    
  # Determine notation to use in output (consistent with what was supplied by the user)
  index.par <- ifelse( xi.notation, 
                       "xi", 
                       "p")
  
  # Remove values of p/xi that are not possible: p <= 1, except p = 0 (the
  # normal distribution). (p = 1 is excluded too: the Poisson likelihood is
  # only defined when every y/phi is an integer, so is not comparable.)
  xi.vec <- xi.vec[ !is.na(xi.vec) ]   # is.na() is also TRUE for NaN
  xi.bad <- (xi.vec <= 1) & (xi.vec != 0)
  if ( any(xi.bad) ) {
    xi.vec <- xi.vec[ !xi.bad ]
    if ( length(xi.vec) == 0 ) {
      stop(paste("No values of", index.par, "are possible: values must be zero, or greater than 1.\n"))
    } else {
      cat("Values of", index.par, "less than or equal to 1 (other than zero) have been removed: such values are not possible.\n")
    }
  }
  

  # With exact zeros in the response, only  p = 0  or  1 < p < 2  are possible
  if ( any( Y == 0 ) ) {
    xi.bad <- (xi.vec >= 2)
    if ( any(xi.bad) ) {
      xi.vec <- xi.vec[ !xi.bad ]
      if ( length(xi.vec) == 0 ) {
        stop(paste("When the response variable contains exact zeros, all values of", index.par, "must be between 1 and 2.\n"))
      } else {
        cat("When the response variable contains exact zeros, all values of", index.par, "must be between 1 and 2; other values have been removed.\n")
      }
    }
  }
  
  if ( do.smooth & (length(xi.vec) < 5) ) {
    warning(paste("Smoothing needs at least five values of", index.par,".") )
    do.smooth <- FALSE
    do.ci     <- FALSE
  }
  if ( (conf.level >= 1) | (conf.level <= 0)  ){
    stop("Confidence level must be between 0 and 1.")
  }
  
  if ( !do.smooth & do.ci ) {
    do.ci <- FALSE
    warning("Confidence intervals only computed if  do.smooth=TRUE\n")
  }
  
  
  if ( add0 & !any(xi.vec == 0) ) {
    xi.vec <- c( 0, xi.vec )
  }

  # Some renaming
  ydata   <- Y
  model.x <- X
  
  # Now, fit models!  Need the Tweedie class
  
  # First define the function to *minimize*;
  # since we want a *maximum* likelihood estimator,
  # define the negative.  (With prior weights, observation  i  has
  # dispersion  phi / weights[i].)
  dtweedie_nlogl <- function(phi, y, mu, power) {
    logl <- sum( dtweedie(y     = y, 
                          mu    = mu, 
                          phi   = phi/weights, 
                          power = power,
                          log   = TRUE ) )
    ans <- -2 * logl
    if ( !is.finite( ans ) ) {
      # The log-likelihood cannot be computed at this phi. Return a huge (but
      # finite) value so optimize() moves away from here. (The mean deviance
      # was used before; it is on a different scale, and so can look like a
      # *good* value of the objective.)
      ans <- 1e300
    }
    # Keep the log-likelihood itself: optimize() returns f() at the minimum
    # (attributes included), so it can be used as L without recomputing it
    attr(ans, "logl") <- logl
    ans
  }
  
  
  
  # The saddlepoint (mean deviance) estimate of phi, allowing for weights
  phi_saddle <- function(mu, power) {
    sum( weights * tweedie_dev(y     = ydata, 
                               mu    = mu, 
                               power = power) ) / length( ydata )
  }
  
  
  # Set up some parameters
  xi.len   <- length(xi.vec)
  L        <- rep(NA_real_, xi.len)
  phi.vec  <- rep(NA_real_, xi.len)
  phi.prev <- NA   # The last good estimate of phi: a start for the next one
  
  # Fitted values from the previous value of p, used to start the next
  # glm.fit(). Neighbouring values of p give similar fits, so IRLS then
  # needs only a few iterations. (NULL means: use the family's usual start.)
  
  mu.start <- NULL

  
  # A progress bar (but not when verbose output is printed, as they clash).
  # on.exit() closes it even if an error occurs.
  show.bar <- (verbose == 0) & (xi.len > 1)
  if ( show.bar ) {
    pBar <- utils::txtProgressBar(min   = 0, 
                                  max   = xi.len,
                                  style = 2)
    on.exit( close(pBar), add = TRUE )
  }
  
  for (i in seq_len(xi.len)) {
    if ( show.bar ) utils::setTxtProgressBar(pBar, i)
    if ( verbose > 0) {
      cat( paste(index.par," = ", xi.vec[i], "\n", sep="") )
    }
    
    # Set things up
    p          <- xi.vec[i]
    L.from.phi <- NULL   # Set below if the search for phi also gives L

    if (verbose == 2) cat("* Fitting initial model:")
    
    # Fit the model with given p
    catch.possible.error <- try(
      fit.model <- stats::glm.fit( x       = model.x, 
                                   y       = ydata, 
                                   weights = weights, 
                                   offset  = offset,
                                   control = control,
                                   mustart = mu.start,
                                   family  = statmod::tweedie(var.power  = p, 
                                                              link.power = link.power)),
      silent = TRUE
    )
    
    
    skip.obs <- inherits( catch.possible.error, "try-error" )

    if( skip.obs ) {
      warning(paste("  Problem near ", index.par, " = ", p,"; this error reported:\n     ",
                    catch.possible.error,
                    " Examine the data and function inputs carefully.") )
      
      mu       <- rep(NA, length(ydata) )
      mu.start <- NULL   # After a failed fit, use the family's default start
    } else {
      mu <- fitted( fit.model )
      # Keep these to start the next fit, if valid for any p (a normal fit,
      # p = 0, may give  mu <= 0  with some links)
      mu.start <- if ( all(mu > 0) ) mu else NULL
    }
    
    ### -- Start:  Estimate phi
    if (verbose >= 1) cat("* Phi estimation, method: ", phi.method)
    
    if( skip.obs ) {
      if (verbose >= 1) cat("; but skipped for this obs\n")
    } else if ( p == 0 ) {
      # Normal distribution: the MLE of phi has a closed form (and is also
      # the mean deviance)
      phi.vec[i] <- sum( weights * (ydata - mu) ^ 2 ) / length(ydata)
      if (verbose >= 1) cat(" (closed form): Done (phi =", phi.vec[i], ")\n")
      
    } else if ( phi.method == "mle" ) {        
        if (verbose >= 1) cat(" (using optimize): ")
        
        # Saddlepoint approx of phi:
        phi.saddle <-  phi_saddle(mu, p)
        
        
        # Start from the previous estimate of phi, if there is one;
        # otherwise, from the saddlepoint approximation
        phi.est   <- if ( is.na(phi.prev) ) phi.saddle else phi.prev
        low.limit <- min( 0.001, phi.saddle / 2)
        
        ans <- try( optimize(f        = dtweedie_nlogl, 
                             maximum  = FALSE, 
                             interval = c(low.limit, 10 * phi.est),
                             power    = p, 
                             mu       = mu, 
                             y        = ydata ),
                    silent = TRUE)
        if ( inherits(ans, "try-error") ) {
          # Treat like a failed fit: record NA, and carry on with the next p
          warning(paste("  Problem estimating phi near ", index.par, " = ", p, "; this error reported:\n     ", ans))
          skip.obs <- TRUE
        } else {
          phi.vec[i] <- ans$minimum
          # The log-likelihood at this phi (from dtweedie()) is L itself
          L.from.phi <- attr(ans$objective, "logl")
          if (verbose >= 1) cat(" Done (phi =", phi.vec[i], ")\n")
    }
  } else{ # phi.method=="saddlepoint")
    if (verbose >= 1) cat(" (using mean deviance/saddlepoint): ")
    phi.vec[i] <- phi_saddle(mu, p)
    if (verbose >= 1) cat(" Done (phi =", phi.vec[i], ")\n")
  }
    
  
  
  phi <- phi.vec[i]
  if ( !is.na(phi) ) phi.prev <- phi
  ### -- Done:  estimate phi
    

  # Now determine log-likelihood at this p.
  
  # With prior weights, the dispersion for each observation is  phi/weights
  
    if (verbose >= 1) cat("* Computing the log-likelihood ")
    
    if (verbose >= 1) cat("(method =", method, "):")
    
    if ( skip.obs ) {
      if (verbose >= 1) cat(" but skipped for this obs\n")
      L[i] <- NA
      
    } else if ( !is.null(L.from.phi) && (method != "saddlepoint") ) {
      # Already computed while finding phi (by dtweedie(), which chooses the
      # fastest accurate method for each observation). Recomputing it with
      # method = "inversion" gives the same value (to about 1e-12), but
      # applies the (slower) inversion to every observation.
      if (verbose >= 1) cat(" (from dtweedie(), while finding phi)")
      L[i] <- L.from.phi    
      } else if ( method == "saddlepoint") {
        L[i] <- dtweedie_logl_saddle(y     = ydata, 
                                     mu    = mu, 
                                     power = p, 
                                     phi   = phi/weights, 
                                     eps   = eps)
      } else if (p == 0) {
        L[i] <- sum( dnorm(x    = ydata, 
                           mean = mu, 
                           sd   = sqrt(phi / weights), 
                           log  = TRUE) )
        
      } else if ( p == 2 ) {
        L[i] <- sum( dgamma(x     = ydata,
                            shape = weights / phi, 
                            rate  = weights / (phi * mu), 
                            log   = TRUE ) )
        
      } else {
        L[i] <- switch(method,
                       "interpolation" = dtweedie_logl( mu    = mu, 
                                                        power = p, 
                                                        phi   = phi / weights, 
                                                        y     = ydata),
                       "series"        = sum( log( dtweedie_series( y     = ydata, 
                                                                    mu    = mu, 
                                                                    power = p, 
                                                                    phi   = phi / weights) ) ),
                       "inversion"     = sum( log( dtweedie_inversion( y     = ydata, 
                                                                       mu    = mu, 
                                                                       power = p, 
                                                                       phi   = phi / weights) ) ) )
        
  }
    
    if (verbose >= 1) {
      cat(" Done: L =", L[i], "\n")
    } 
  }
  

  ### Only finite values of L are used for smoothing and finding the maximum
  # (is.finite() is FALSE for NA, NaN and +/-Inf)
  ok          <- is.finite(L)
  xi.vec.fix  <- xi.vec[ok]
  L.fix       <- L[ok]
  phi.vec.fix <- phi.vec[ok]
  
  if ( !any(ok) ) {
    stop("No valid values of the log-likelihood were computed.\n",
         "  Consider trying other values of ", index.par, ", or another value for  method.")
  }
  if ( any(!ok) ) {
    cat(" WARNING: Some values of L are infinite or NA, and are ignored\n")
  }
  if ( do.smooth & (length(L.fix) < 5) ) {
    warning(paste("Smoothing needs at least five (valid) values of", index.par, "; no smoothing used.") )
    do.smooth <- FALSE
    do.ci     <- FALSE
  }


  ### Smoothing and plotting
  # x and y are returned (to recreate the plot): the smoothed profile when
  # smoothing, otherwise the computed values
  
  # With add0, the spline would join p = 0 to the values p > 1, across the 
  # impossible values 0 < p < 1. So the MLE and CI use a spline through the
  # values p > 1 only (ss.pos), and p = 0 is then considered on its own.
  xi.pos <- (xi.vec.fix > 0)
  
  if ( do.smooth ) {
    if (verbose >= 1) cat(" --- \n* Smoothing: ")
    # An interpolating spline through the points (it does not smooth them)
    ss     <- stats::splinefun( xi.vec.fix, 
                                L.fix )
    ss.pos <- stats::splinefun( xi.vec.fix[xi.pos],
                                L.fix[xi.pos] )
    
    x <- seq(min(xi.vec.fix), 
             max(xi.vec.fix), 
             length = 50 )
    y <- ss( x )
  } else {
    x <- xi.vec
    y <- L
  }

  if ( do.plot ) {
    # Plot only the finite values
    x.plot <- if ( do.smooth ) x else xi.vec.fix
    y.plot <- if ( do.smooth ) y else L.fix
  
        do.call(plot,
                c(list(
                    x    = range(xi.vec),
                    y    = range(y.plot),
                    type = "n",
                    las  = 1,
                    xlab = ifelse(xi.notation, 
                                  expression(paste( xi, " index")),
                                  expression(paste( italic(p), " index")) ),
                    ylab = expression(italic(L))
                ),
                plot_args))

        do.call(lines,
                c(list(
                    x = x.plot,
                    y = y.plot),
                    line_args))
        
        rug( xi.vec )
        if (do.points) {
          do.call(points,
                  c(list(
                      x = xi.vec,
                      y = L),
                     point_args))
        }
        if (add0 & do.smooth) lines(x[x < 1],
                                    y[x < 1],
                        col = "gray",
                        lwd = 2)
      }
  
  if (verbose >= 2) cat(" Done\n")
  

  ### Maximum likelihood estimates
  at.end <- FALSE   # Is the maximum at the end of the values of p?

  if ( do.smooth ){
    if (verbose >= 2) cat(" Finding maximum likelihood estimates:  ")
    
    # Find the best point on a grid, then refine with optimize() between
    # its neighbours (the grid alone is accurate only to range/49)
    x.grid <- seq(min(xi.vec.fix[xi.pos]), 
                  max(xi.vec.fix[xi.pos]), 
                  length = 50)
    y.grid <- ss.pos(x.grid)
    best   <- which.max(y.grid)
    
    if ( (best == 1) | (best == length(x.grid)) ) {
      at.end <- TRUE
      xi.max <- x.grid[best]
      L.max  <- y.grid[best]
    } else {
      opt    <- optimize(ss.pos, 
                         interval = x.grid[c(best - 1, best + 1)], 
                         maximum  = TRUE)
      xi.max <- opt$maximum
      L.max  <- opt$objective
    }
    
    # With add0, p = 0 (the normal distribution) may be better still
    if ( any(xi.vec.fix == 0) ) {
      L.0 <- L.fix[xi.vec.fix == 0]
      if ( L.0 > L.max ) {
        xi.max <- 0
        L.max  <- L.0
        at.end <- TRUE
      }
    }
    
    # Now find the mle of  phi  at this very value of  p
    if ( at.end ) {
      # At one of the computed values: use its phi
      phi.max <- phi.vec.fix[ xi.vec.fix == xi.max ]
      if (xi.max != 0) warning("True maximum possibly not detected.")
      
    } else {
      fit.max <- stats::glm.fit( y       = ydata, 
                                 x       = model.x, 
                                 weights = weights, 
                                 offset  = offset,
                                 control = control,
                                 family  = statmod::tweedie(var.power  = xi.max, 
                                                            link.power = link.power))
      mu <- fitted( fit.max )
      
      if ( phi.method == "saddlepoint"){
        phi.max <- phi_saddle(mu, xi.max)
      } else {
        # Bound phi using the estimates at the computed values of p
        # (only p > 1: phi at p = 0 is on a different scale)
        phi.lo <- 0.5 * min( phi.vec.fix[xi.pos] )
        phi.hi <-   2 * max( phi.vec.fix[xi.pos] )
        
        phi.max <- optimize( f        = dtweedie_nlogl, 
                             maximum  = FALSE, 
                             interval = c(phi.lo, phi.hi ), 
                             power    = xi.max, 
                             mu       = mu,
                             y        = ydata)$minimum
      }
    }
    if (verbose >= 2) cat(" Done\n")
  } else {
    # Find maximum from computed data:  which.max()  takes the first if
    # there are ties
    if (verbose >= 2) cat(" Finding maximum likelihood estimates:  ")
    best    <- which.max(L.fix)
    xi.max  <- xi.vec.fix[best]
    L.max   <- L.fix[best]
    phi.max <- phi.vec.fix[best]
    
    if (verbose >= 2) cat(" Done\n")
  }
  
  # Now report
  if ( verbose >= 2 ) {
    cat( "ML Estimates:  ", index.par, "=",xi.max, " with phi=", phi.max, " giving L=", L.max, "\n")
    cat(" ---\n")
  }
  
  # Now find approximate, nominal 95% CI info
  ht <- L.max - ( stats::qchisq(conf.level, 1) / 2 )
  ci <- c(NA, NA )
  
  
  if ( do.ci ) {
    if (verbose == 2) cat("* Finding confidence interval:")

    if ( do.plot ) {
      graphics::abline(h   = ht, 
                       lty = 2)
      graphics::title( sub = paste("(", 100*conf.level, "% confidence interval)",
                       sep = "") )
    }
    
    # Each limit is where the spline crosses the height  ht: solve
    #   ss.pos(p) = ht  
    # using uniroot(), on each side of the maximum
    ci.fun <- function(z) ss.pos(z) - ht
    p.lo   <- min( x.grid )
    p.hi   <- max( x.grid )
    
    if ( xi.max < p.lo ) {
      warning("Confidence interval not found: the MLE is p = 0.\n")
    } else {
      # --- Left side ---
      if ( ci.fun(p.lo) >= 0 ) {
        warning("Confidence interval cannot be found: insufficient data to find left CI.\n")
      } else {
        ci[1] <- stats::uniroot(ci.fun, 
                                interval = c(p.lo, xi.max))$root
      }
    # --- Right side ---
      if ( ci.fun(p.hi) >= 0 ) {
        warning("Confidence interval cannot be found: insufficient data to find right CI.\n")
    } else {
      
      ci[2] <- stats::uniroot(ci.fun, 
                              interval = c(xi.max, p.hi))$root
    }
      }
    if (verbose==2) cat(" Done\n")
  }
  
    
    # Build the output (in the notation the user supplied)
    
      out <- list( y          = y, 
                   x          = x, 
                   ht         = ht, 
                   L          = L, 
                   xi         = xi.vec,
                   xi.max     = xi.max, 
                   L.max      = L.max, 
                   phi        = phi.vec, 
                   phi.max    = phi.max, 
                   ci         = ci,
                   method     = method, 
                   phi.method = phi.method)
      
      if ( fit.glm ) {
        out$glm.obj <- stats::glm.fit( x       = model.x, 
                                       y       = ydata, 
                                       weights = weights, 
                                       offset  = offset,
                                       control = control,
                                       family  = statmod::tweedie(var.power  = xi.max, 
                                                                  link.power = link.power) )
      }
      
      
      if ( !xi.notation ) {
        names(out)[ names(out) == "xi" ]     <- "p"
        names(out)[ names(out) == "xi.max" ] <- "p.max"
      }
  
  if ( verbose ) cat("\n")
  
  invisible( out )
  
}




#' @rdname tweedie_profile
#' @export
tweedie.profile <- function(formula, 
                            p.vec      = NULL, 
                            xi.vec     = NULL, 
                            link.power = 0, 
                            data, 
                            weights    = 1, 
                            offset     = 0, 
                            fit.glm    = FALSE, 
                            do.smooth  = TRUE, 
                            do.plot    = FALSE, 
                            do.ci      = do.smooth, 
                            eps        = 1/6,
                            control    = list( epsilon = 1e-09, 
                                               maxit   = stats::glm.control()$maxit, 
                                              trace    = glm.control()$trace ),
                            do.points  = do.plot, 
                            method     = "inversion",
                            conf.level = 0.95, 
                            phi.method = ifelse(method == "saddlepoint", "saddlepoint", "mle"), 
                            verbose    = FALSE, 
                            add0       = FALSE,
                            plot_args  = list(),
                            point_args = list(),
                            line_args  = list()) {  
  # 1. Signal the deprecation
  lifecycle::deprecate_warn(
    when = "3.0.5", 
    what = "tweedie.profile()", 
    with = "tweedie_profile()"
  )
  
  # 2. Capture the call as the user wrote it
  cl <- match.call()
  
  # 3. Change the name of the function to be called
  cl[[1]] <- quote(tweedie::tweedie_profile)
  
  # 4. Execute it in the environment where tweedie.profile was called
  eval(cl, parent.frame())
}


