\name{networkScreeningGS}
\alias{networkScreeningGS}
\title{ Network gene screening with an external gene significance measure }
\description{
  This function blends standard and network approaches to selecting genes (or variables in general) with
high gene significance
}
\usage{
networkScreeningGS(
  datExpr,
  datME,
  GS,
  oddPower = 3,
  blockSize = 1000,
  minimumSampleSize = ..minNSamples, 
  addGS = TRUE)
}
\arguments{
  \item{datExpr}{ data frame of expression data }
  \item{datME}{ data frame of module eigengenes }
  \item{GS}{ numeric vector of gene significances }
  \item{oddPower}{ odd integer used as a power to raise module memberships and significances }
  \item{blockSize}{ block size to use for calculations with large data sets }
  \item{minimumSampleSize}{ minimum acceptable number of samples. Defaults to the default minimum number of
samples used throughout the WGCNA package, currently 4.}
  \item{addGS}{ logical: should gene significances be added to the screening statistics?}
}
\details{
This function should be considered experimental. It takes into account both the "standard" and the network
measures of gene importance for the trait.
}
\value{
  \item{GS.Weighted }{weighted gene significance }
  \item{GS }{copy of the input gene significances (only if \code{addGS=TRUE})}
}
\author{ Steve Horvath }
\seealso{\code{\link{networkScreening}}, \code{\link{automaticNetworkScreeningGS}}}
\keyword{ misc}
