\name{robustscale}
\alias{robustscale}
\title{Robust Estimation of Median (center) and MAD (scale)}
\description{
    Smoothing spline estimate of median and mean absolute deviation (MAD).
}
\usage{
robustscale(y, x, nslice=400, corcenter=TRUE, decrease=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{response}
  \item{x}{predictor}
  \item{nslice}{number of slices (should be "large")}
  \item{corcenter}{correct for center}
  \item{decrease}{force MAD to decrease with \code{x}}
}
\details{
    This divides data into roughly many \code{nslice} slices and computes
    median and mean absolute deviation (\code{mad}) for each slice. These are
    then smoothed using \code{smooth.spline}.
}
\value{
    Data frame containing significant genes with the following information:
  \item{center}{estimate of center median}
  \item{scale}{MAD estimate of scale}
  \item{x}{ordered \code{x} values for plotting}
  \item{y}{\code{y} sorted by \code{x}}
}
\author{Yi Lin}
\seealso{\code{\link{mad}}, \code{\link{smooth.spline}}}
\examples{
\dontrun{
robustscale(y,x)
}
}
\keyword{robust}
\keyword{smooth}
