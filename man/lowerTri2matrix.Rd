\name{lowerTri2matrix}
\alias{lowerTri2matrix}
\title{
Reconstruct a symmetric matrix from a distance (lower-triangular) representation
}
\description{
Assuming the input vector contains a vectorized form of the distance representation of a symmetric matrix,
this function creates the corresponding matrix. This is useful when re-forming symmetric matrices that have
been vectorized to save storage space.
}
\usage{
lowerTri2matrix(x, diag = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a numeric vector 
}
  \item{diag}{
 value to be put on the diagonal. Recycled if necessary.
}
}
\details{
The function assumes that \code{x} contains the vectorized form of the distance representation of a
symmetric matrix. In particular, \code{x} must have a length that can be expressed as n*(n-1)/2, with n an
integer. The result of the function is then an n times n matrix. 
}
\value{
  A symmetric matrix whose lower triangle is given by \code{x}.
}
\author{
Peter Langfelder
}
\examples{
  # Create a symmetric matrix
  m = matrix(c(1:16), 4,4)
  mat = (m + t(m));
  diag(mat) = 0;

  # Print the matrix
  mat

  # Take the lower triangle and vectorize it (in two ways)
  x1 = mat[lower.tri(mat)]
  x2 = as.vector(as.dist(mat))

  all.equal(x1, x2) # The vectors are equal

  # Turn the vectors back into matrices
  new.mat = lowerTri2matrix(x1, diag = 0);

  # Did we get back the same matrix?

  all.equal(mat, new.mat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
