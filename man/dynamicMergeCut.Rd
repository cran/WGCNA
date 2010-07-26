\name{dynamicMergeCut}
\alias{dynamicMergeCut}
\title{ Threshold for module merging }
\description{
  Calculate a suitable threshold for module merging based on the number of samples and a desired Z
quantile.
}
\usage{
dynamicMergeCut(n, mergeCor = 0.9, Zquantile = 2.35)
}
\arguments{
  \item{n}{number of samples }
  \item{mergeCor}{ theoretical correlation threshold for module merging  }
  \item{Zquantile}{ Z quantile for module merging }
}
\details{
  This function calculates the threshold for module merging. The threshold is calculated as the lower
boundary of the interval around the theoretical correlation \code{mergeCor} whose width is given by the
Z value \code{Zquantile}.
}
\value{
  The correlation threshold for module merging; a single number.
}
\author{ Steve Horvath }
\seealso{ \code{\link{moduleEigengenes}}, \code{\link{mergeCloseModules}} }
\examples{
  dynamicMergeCut(20)
  dynamicMergeCut(50)
  dynamicMergeCut(100)
}
\keyword{ misc }% __ONLY ONE__ keyword per line
