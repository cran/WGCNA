\name{removePrincipalComponents}
\alias{removePrincipalComponents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove leading principal components from data
}
\description{
This function calculates a fixed number of the first principal components of the given data and returns the
residuals of a linear regression of each column on the principal components. 
}
\usage{
removePrincipalComponents(x, n)
}
\arguments{
  \item{x}{
Input data, a numeric matrix. All entries must be non-missing and finite.
}
  \item{n}{
Number of principal components to remove. This must be smaller than the smaller of the number of rows and
columns in \code{x}.
}
}
\value{
A matrix of residuals of the same dimensions as \code{x}.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{svd}} for singular value decomposition, 
\code{\link{lm}} for linear regression
}
\keyword{misc}
