\name{oddsplot}
\alias{oddsplot}
\title{Odds Plot for Differential Microarray Expression}
\description{
    The function plots contours for the odds that points on microarray show
    differential expression between two conditions (e.g. Cy3 and Cy5 dye
    channels on the same microarray).
}
\usage{
oddsplot(x, y, theta, by.level = 10, rotate = FALSE, offset =
         0, main = "", xlab = xlabs, ylab = ylabs, col = NULL,
         cex = c(0.25, 0.75), shrink = FALSE, lims =
         range(c(x, y)))
}
\arguments{
  \item{x}{first condition expression levels}
  \item{y}{second condition expression levels}
  \item{theta}{four parameters from \code{em.ggb}}
  \item{by.level}{odds plot contours increase by this level}
  \item{rotate}{rotate to average versus ratio if TRUE, otherwise plot
    conditions against each other}
  \item{offset}{offset for \code{log} transform}
  \item{main}{main title for plot}
  \item{xlab}{horizontal axis label (default if \code{Cy3} if
      \code{rotate} is FALSE, \code{Average Intensity} otherwise}
  \item{ylab}{vertical axis label (default if \code{Cy5} if
      \code{rotate} is FALSE, \code{Cy3 / Cy5} otherwise}
  \item{col}{color of points (if NULL, use black for non-changing
      points, blue for changing points)}
  \item{cex}{character expansion (use \code{rep(.25,2)} to have all points the
      same size)}
  \item{shrink}{use shrinkage on expression levels if TRUE (default is FALSE)}
  \item{lims}{limits for plot area}
}
\details{
Fit Gamma/Gamma/Bernoulli model (equal marginal distributions)
The model has spot intensities x ~ Gamma(a,b); y ~ Gamma(a,c).
The shape parameters b and c are ~ Gamma(a0,nu).
With probability p, b = c; otherwise b != c. All spots are assumed to be
independent.}
\value{
    Log odds for all points in original order.
}
\references{MA Newton, CM Kendziorski, CS Richmond, FR Blattner and KW
    Tsui (2000) ``On differential variability of expression ratios:
    improving statistical inference about gene expression changes from
    microarray data,''
    \emph{J Computational Biology 00}: 000-000.
}
\author{Michael Newton}
\seealso{\code{\link{em.ggb}}}
\examples{
\dontrun{
oddsplot( x, y )
}}
\keyword{hplot}
\keyword{models}
