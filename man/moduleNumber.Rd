\name{moduleNumber}
\alias{moduleNumber}
\title{Fixed-height cut of a dendrogram.}
\description{
  Detects branches of on the input dendrogram by performing a fixed-height cut.
}
\usage{
moduleNumber(dendro, cutHeight = 0.9, minSize = 50)
}
\arguments{
  \item{dendro}{a hierarchical clustering dendorgram such as one returned by \code{hclust}.  }
  \item{cutHeight}{Maximum joining heights that will be considered.  }
  \item{minSize}{Minimum cluster size.  }
}
\details{
  All contiguous branches below the height \code{cutHeight} that contain at least \code{minSize} objects
are assigned unique positive numerical labels; all unassigned objects are assigned label 0.
}
\value{
  A vector of numerical labels giving the assigment of each object.
}
\note{The numerical labels may not be sequential. See \code{\link{normalizeLabels}} for a way to put the
labels into a standard order.}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{
  \code{\link{hclust}}, \code{\link{cutree}}, \code{\link{normalizeLabels}}
}
\keyword{cluster}
