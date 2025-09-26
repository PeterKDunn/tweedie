\name{tweedie integrand}
\alias{tweedie_Igrand}
\title{Tweedie Distributions: density integrand}
\description{The integrand of the Fourier integral for evaluating Tweedie density function}
\usage{
  tweedie_Igrand(phi, power, y, mu, t, doPlot, type)
}
\arguments{
  \item{phi}{the dispersion}
  \item{power}{the value of \eqn{p}{power} such that the variance is \eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}}
  \item{y}{the value at which to evaluate the density}
  \item{mu}{the value of the mean; when \code{type} is \code{"PDF"}, the density is always evaluatred at \code{mu = 1}, so the value of \code{mu} given makes no difference}
  \item{t}{the range over which to evaluate the integrand}
  \item{doPlot}{if \code{TRUE} (the default), a plot of the integrand is produced}
  \item{type}{if \code{PDF}, the integrand for the density is produced; otherwise, the integrand for teh DF is produced}
}
\value{
  \code{ig}, a vector, of the same length as \code{y}, containing the evaluation of the integrand for the Fourier inversion to compute the density;
  \code{k.real}, the real component of the integrand;
  \code{k.imag}, the real component of the integrand. 
  For \eqn{p>2}{power > 2}, the integrand is \eqn{\exp{\Re k}\cos{\Im k}}{exp( k.real ) * cos( k.imag)} for the density.
  When \eqn{1<p<2}{1 < power < 2}, there needs to be adjustment for working with just the continuous part; so when \eqn{1<p<2}{1 < power < 2}, the integrand is
  \eqn{\exp{\Re k}\cos{\Im k} - \exp( -\lambda_1)\cos(ty)}{exp( k.real ) * cos( k.imag) - exp( lambda1 ) * cos(ty)},
  where \eqn{\lambda_1 = 1/[\phi(2 - p)]}{lambda1 = 1/[ phi*(2 - power) ]}.

  The distribution function emerges along similar lines, but has a \eqn{1/t}{1/t} term in the integrand, so converges faster than the density, in general.
  
  Note that \code{ig} is just the integrand; it does not include the factors in front of the integrand.
  In all cases, we are assuming \eqn{\mu = 1}{mu = 1}.
}
\details{
  The integrand appears in the Fourier inversion of the Tweedie densities.
  For details, see \code{?dtweedie.inversion} and the references therein.
}
\author{Peter Dunn (\email{dunn@usq.edu.au})}
\references{
Dunn, Peter K and Smyth, Gordon K (To appear).
Series evaluation of Tweedie exponential dispersion model densities
\emph{Statistics and Computing}.

Dunn, Peter K and Smyth, Gordon K (2001).
Tweedie Family Densities: Methods of Evaluation
\emph{Proceedings of the 16th International Workshop on Statistical Modelling},
Odense, Denmark, 2--6 July

Jorgensen, B. (1987).
Exponential dispersion models
\emph{Journal of the Royal Statistical Society}, B,
\bold{49}, 127--162.

Jorgensen, B. (1997).
\emph{Theory of Dispersion Models},
Chapman and Hall, London.

Sidi, Avram (1982).
The numerical evaluation of very oscillatory infinite integrals by
extrapolation
\emph{Mathematics of Computation}
\bold{38}(158), 517--529.

Sidi, Avram (1988).  
A user-friendly extrapolation method for
oscillatory infinite integrals
\emph{Mathematics of Computation}
\bold{51}(183), 249--266.

Tweedie, M. C. K. (1984).
An index which distinguishes between some important exponential families
\emph{Statistics: Applications and New Directions.
Proceedings of the Indian Statistical Institute Golden Jubilee International Conference}
(Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
}
\seealso{\code{dtweedie}, \code{dtweedie.inversion} in the \code{tweedie} package }

\keyword{models}

