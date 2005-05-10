\name{pickgene}
\alias{pickgene}
\title{Plot and Pick Genes based on Differential Expression}
\description{
    The function picks plots the average intensity versus linear
    contrasts (currently linear, quadratic up to cubic) across
    experimental conditions. Critical line
    is determine according to Bonferroni-like multiple comparisons,
    allowing SD to vary with intensity.
}
\usage{
pickgene(data, geneID = 1:nrow(data), overalllevel = 0.05,
         npickgene = -1, marginal = FALSE, rankbased = TRUE,
         allrank = FALSE, meanrank = FALSE, offset = 0,
         modelmatrix = model.pickgene(faclevel, facnames,
         contrasts.fac, collapse, show, renorm), faclevel =
         ncol(data), facnames =
         letters[seq(length(faclevel))], contrasts.fac =
         "contr.poly", show = NULL, main = "", renorm = 1,
         drop.negative = FALSE, plotit = npickgene < 1, mfrow
         = c(nr, nc), mfcol = NULL, ylab = paste(shownames,
         "Trend"), ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{data matrix}
  \item{geneID}{gene identifier (default \code{1:nrow(x)})}
  \item{overalllevel}{overall significance level (default \code{0.05})}
  \item{npickgene}{number of genes to pick (default \code{-1} allows
      automatic selection)}
  \item{marginal}{additive model if TRUE, include interactions if FALSE}
  \item{rankbased}{use ranks if TRUE, log tranform if FALSE}
  \item{allrank}{rank all chips together if true, otherwise rank separately}
  \item{meanrank}{show mean abundance as rank if TRUE}
  \item{offset}{offset for log transform}
  \item{modelmatrix}{model matrix with first row all 1's and other rows
      corresponding to design contrasts; automatically created by call
      to \code{model.pickgene} if omitted}
  \item{faclevel}{number of factor levels for each factor}
  \item{facnames}{factor names}
  \item{contrasts.fac}{type of contrasts}
  \item{show}{vector of contrast numbers to show (default is all)}
  \item{main}{vector of main titles for plots (default is none)}
  \item{renorm}{vector to renormalize contrasts (e.g. use \code{sqrt(2)}
      to turn two-condition contrast into fold change)}
  \item{drop.negative}{drop negative values in log transform}
  \item{plotit}{plot if TRUE}
  \item{mfrow}{\code{par()} plot arrangement by rows (default up to 6
      per page; set to NULL to not change)}
  \item{mfcol}{\code{par()} plot arrangement by columns (default is NULL)}
  \item{ylab}{vertical axis labels}
  \item{...}{parameters for \code{robustscale}}
}
\details{
Infer genes that differentially express across conditions using a robust
data-driven method. Adjusted gene expression levels \code{A} are
replaced by \code{qnorm(rank(A))}, followed by \code{robustscale}
estimation of center and spread. Then Bonferroni-style gene by gene
tests are performed and displayed graphically.
}
\value{
    Data frame containing significant genes with the following information:
    \item{pick}{data frame with picked genes}
    \item{score}{data frame with center and spread for plotting}
    Each of these is a list with elements for each contrast.
    The \code{pick} data frame elements have the following information:
  \item{probe}{gene identifier}
  \item{average}{average gene intensity}
  \item{fold1}{positive fold change}
  \item{fold2}{negative fold change}
  \item{pvalue}{Bonferroni-corrected p-value}
  The \code{score} data frame elements have the following:
  \item{x}{mean expression level (antilog scale)}
  \item{y}{contrast (antilog scale)}
  \item{center}{center for contrast}
  \item{scale}{scale (spread) for contrast}
  \item{lower}{lower test limit}
  \item{upper}{upper test limit}
}
\references{Y Lin, BS Yandell and ST Nadler (2000)
    ``Robust Data-Driven Inference for Gene Expression Microarray
    Experiments,'' Technical Report, Department of Statistics, UW-Madison.
}
\author{Yi Lin and Brian Yandell}
\seealso{\code{\link{pickgene}}}
\examples{
\dontrun{
pickgene( data )
}
}
\keyword{hplot}
\keyword{models}

