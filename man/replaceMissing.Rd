\name{replaceMissing}
\alias{replaceMissing}
\title{
Replace missing values with a constant.
}
\description{
A convenience function for replacing missing values with a (non-missing) constant.
}
\usage{
replaceMissing(x, replaceWith)
}
\arguments{
  \item{x}{
An atomic vector or array.
}
  \item{replaceWith}{
Value to replace missing entries in \code{x}. The default is \code{FALSE} for logical vectors, 
0 for numeric vectors, and empty string "" for character vectors.
}
}
\value{
\code{x} with missing data replaced.
}
\author{
Peter Langfelder
}
\examples{
logVec = c(TRUE, FALSE, NA, TRUE);
replaceMissing(logVec)

numVec = c(1,2,3,4,NA,2)
replaceMissing(numVec)

}
\keyword{misc}
