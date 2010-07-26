\name{checkAdjMat}
\alias{checkAdjMat}
\title{ Check adjacency matrix }
\description{
  Checks a given matrix for properties that an adjacency matrix must satisfy.
}
\usage{
checkAdjMat(adjMat, min = 0, max = 1)
}
\arguments{
  \item{adjMat}{ matrix to be checked }
  \item{min}{minimum allowed value for entries of adjMat}
  \item{max}{maximum allowed value for entries of adjMat}
}
\details{
  The function checks whether the given matrix really is a 2-dimensional numeric matrix, whether it is
square, symmetric, and all finite entries are between \code{min} and \code{max}. 
If any of the conditions is not met, the
function issues an error. 
}
\value{
  None. The function returns normally if all conditions are met.
}
\author{Peter Langfelder}
\seealso{\code{\link{adjacency}}}
\keyword{misc }
