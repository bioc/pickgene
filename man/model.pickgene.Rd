\name{model.pickgene}
\alias{model.pickgene}
\title{Create Model Matrix for Orthogonal Contrasts}
\description{
    The function created a model matrix of orthogonal contrasts to be
    used by pickgene.
}
\usage{
model.pickgene(faclevel, facnames = letters[seq(length(faclevel))],
               contrasts.fac = "contr.poly", collapse = "+", show =
               NULL, renorm = 1, modelexpr = formula(paste("~",
               paste(facnames, collapse = collapse))),
               contrasts.list = contr.list)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{faclevel}{vector with number of levels for each factor}
  \item{facnames}{vector of factor names (default = "a","b",...)}
  \item{contrasts.fac}{vector of contrast types}
  \item{collapse}{"+" for additive model, "*" for full model with interactions}
  \item{show}{vector of contrast numbers to show (default is all)}
  \item{renorm}{vector to renormalize contrasts (e.g. use \code{sqrt(2)}
      to turn two-condition contrast into fold change)}
  \item{modelexpr}{model formula}
  \item{contrasts.list}{list of contrasts indexed by \index{facnames}}
}
\details{
Creates a model matrix data frame with first column having all 1's and
other columns having contrasts.
}
\value{
Result of call to model.matrix
}
\author{Brian Yandell}
\seealso{\code{\link{model.matrix}}}
\examples{
model.pickgene( c(2,3), c("sex","genotype"))
}
\keyword{utilities}
