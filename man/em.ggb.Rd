\name{em.ggb}
\alias{em.ggb}
\title{EM calculation for Gamma-Gamma-Bernoulli Model}
\description{
    The function plots contours for the odds that points on microarray show
    differential expression between two conditions (e.g. Cy3 and Cy5 dye
    channels on the same microarray).
}
\usage{
em.ggb(x, y, theta, start = c(2,1.2,2.7), pprior = 2,
  printit = FALSE, tol = 1e-9, offset = 0 )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{first condition expression levels}
  \item{y}{second condition expression levels}
  \item{theta}{four parameters \code{a,a0,nu,p}}
  \item{start}{starting estimates for theta}
  \item{pprior}{Beta hyperparameter for prob \code{p} of differential
      expression} 
  \item{printit}{print iterations if TRUE}
  \item{tol}{parameter tolerance for convergence}
  \item{offset}{offset added to xx and yy before taking log (can help
      with negative adjusted values)}
}
\details{
Fit Gamma/Gamma/Bernoulli model (equal marginal distributions)
The model has spot intensities x ~ Gamma(a,b); y ~ Gamma(a,c).
The shape parameters b and c are ~ Gamma(a0,nu).
With probability p, b = c; otherwise b != c. All spots are assumed to be
independent.}
\value{
    Four parameter vector \code{theta} after convergence.
}
\references{MA Newton, CM Kendziorski, CS Richmond, FR Blattner and KW
    Tsui (2000) ``On differential variability of expression ratios:
    improving statistical inference about gene expression changes from
    microarray data,''
    \emph{J Computational Biology 00}: 000-000.
}
\author{Michael Newton}
\seealso{\code{\link{oddsplot}}, \code{\link{Ecoli.oddsplot}}}
\examples{
\dontrun{
em.ggb( x, y )
}}
\keyword{models}
