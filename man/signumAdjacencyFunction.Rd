\name{signumAdjacencyFunction}
\alias{signumAdjacencyFunction}

\title{ Hard-thresholding adjacency function }
\description{
  This function transforms correlations or other measures of similarity into an unweighted network
adjacency.
}
\usage{
signumAdjacencyFunction(corMat, threshold)
}
\arguments{
  \item{corMat}{ a matrix of correlations or other measures of similarity. }
  \item{threshold}{ threshold for connecting nodes: all nodes whose \code{corMat} is above the threshold
will be connected in the resulting network. }
}
\value{
  An unweighted adjacency matrix of the same dimensions as the input \code{corMat}.
}
\references{
      Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 
}

\author{ Steve Horvath }
\seealso{ \code{\link{adjacency}} for soft-thresholding and creating weighted networks. }
\keyword{ misc }% __ONLY ONE__ keyword per line
