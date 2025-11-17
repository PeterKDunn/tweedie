#' Tweedie internal function
#'
#' @name Tweedie internals
#' @aliases dtweedie.dlogfdphi dtweedie.logl dtweedie.logl.saddle dtweedie.logv.bigp dtweedie.logw.smallp dtweedie.interp dtweedie.jw.smallp dtweedie.kv.bigp dtweedie.series.bigp dtweedie.series.smallp stored.grids check_Inputs sort_Notation special_Cases_CDF dtweedie_Fortran dtweedie_Inversion_Report dtweedie_Inversion_Threemethods dtweedie_Series_Report ptweedie_Inversion_Report ptweedie_Series_Report special_Cases_CDF
#' @title Tweedie internal function
#' @description Internal tweedie functions. These are not to be called by the user.
#'
#' @usage
#' dtweedie.dlogfdphi(y, mu, phi, power)
#' dtweedie.logl(phi, y, mu, power)
#' dtweedie.logl.saddle( phi, power, y, mu, eps=0)
#' dtweedie.logv.bigp( y, phi, power)
#' dtweedie.logw.smallp(y, phi, power)
#' dtweedie.interp(grid, nx, np, xix.lo, xix.hi,p.lo, p.hi, power, xix)
#' dtweedie.jw.smallp(y, phi, power )
#' dtweedie.kv.bigp(y, phi, power)
#' dtweedie.series.bigp(power, y, mu, phi)
#' dtweedie.series.smallp(power, y, mu, phi)
#' stored.grids(power)
#' 
#'
#'
#' @param y the vector of responses
#' @param y the vector of responses
#' @param power the value of \eqn{p}{power} such that the variance is \eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}
#' @param mu the mean
#' @param phi the dispersion
#' @param grid the interpolation grid necessary for the given value of \eqn{p}{power}
#' @param nx the number of interpolation points in the \eqn{\xi}{xi} dimension
#' @param np the number of interpolation points in the \eqn{p}{power} dimension
#' @param xix.lo the lower value of the transformed \eqn{\xi}{xi} value used in the interpolation grid. (Note that the value of \eqn{\xi}{xi} is from \eqn{0} to \eqn{\infty}{infty}, and is transformed such that it is on the range \eqn{0} to \eqn{1}.)
#' @param xix.hi the higher value of the transformed \eqn{\xi}{xi} value used in the interpolation grid.
#' @param p.lo the lower value of the \eqn{p} value used in the interpolation grid.
#' @param p.hi the higher value of the \eqn{p} value used in the interpolation grid.
#' @param xix the value of the transformed \eqn{\xi}{xi} at which a value is sought.
#' @param eps the offset in computing the variance function in the saddlepoint approximation. The default is \code{eps=1/6} (as suggested by Nelder and Pregibon, 1987).
#' @param p the Tweedie index parameter
#'
#' @author Peter Dunn (\email{pdunn2@usc.edu.au})
#' @references
#' 	Nelder, J. A. and Pregibon, D. (1987).
#' 	An extended quasi-likelihood function
#' 	\emph{Biometrika},
#' 	\bold{74}(2), 221--232.
#' 	\doi{10.1093/biomet/74.2.221}
#' @keywords models
NULL