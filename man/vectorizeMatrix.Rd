\name{vectorizeMatrix}
\alias{vectorizeMatrix}
\title{ Turn a matrix into a vector of non-redundant components }
\description{
  A convenient function to turn a matrix into a vector of non-redundant components. If the matrix is
non-symmetric, returns a vector containing all entries of the matrix. If the matrix is symmetric, only
returns the upper triangle and optionally the diagonal.
}
\usage{
vectorizeMatrix(M, diag = FALSE)
}
\arguments{
  \item{M}{ the matrix or data frame to be vectorized. }
  \item{diag}{ logical: should the diagonal be included in the output? }
}
\value{
  A vector containing the non-redundant entries of the input matrix.
}
\author{ Steve Horvath }
\keyword{ misc }% __ONLY ONE__ keyword per line
