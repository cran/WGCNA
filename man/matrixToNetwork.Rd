\name{matrixToNetwork}
\alias{matrixToNetwork}
\title{
  Construct a network from a matrix
}
\description{
Constructs a network
}
\usage{
matrixToNetwork(
    mat, 
    symmetrizeMethod = c("average", "min", "max"), 
    signed = TRUE, 
    min = NULL, max = NULL, 
    power = 12,
    diagEntry = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{matrix to be turned into a network. Must be square.
}
  \item{symmetrizeMethod}{
  method for symmetrizing the matrix. The method will be applied to each component of mat and its transpose.
}
  \item{signed}{
   logical: should the resulting network be signed? Unsigned networks are constructed from \code{abs(mat)}.
}
  \item{min}{
  minimum allowed value for \code{mat}. If \code{NULL}, the actual attained minimum of \code{mat} will be used.
Missing data are ignored. Values below \code{min} are truncated to \code{min}.
}
  \item{max}{
  maximum allowed value for \code{mat}. If \code{NULL}, the actual attained maximum of \code{mat} will be used.
Missing data are ignored. Values below \code{max} are truncated to \code{max}.
}
  \item{power}{
  the soft-thresholding power.
}
  \item{diagEntry}{
  the value of the entries on the diagonal in the result. This is usally 1 but some applications may require
a zero (or even NA) diagonal.
}
}
\details{
  If \code{signed} is \code{FALSE}, the matrix \code{mat} is first converted to its absolute value.

  This function then symmetrizes the matrix using the \code{symmetrizeMethod} component-wise on \code{mat}
and \code{t(mat)} (i.e., the transpose of \code{mat}). 

  In the next step, the symmetrized matrix is linearly scaled to the interval [0,1] using either \code{min}
 and \code{max} (each either supplied or determined from the matrix). Values outside of the [min, max] range
are truncated to \code{min} or \code{max}. 

  Lastly, the adjacency is calculated by rasing the matrix to \code{power}. 
The diagonal of the result is set to
\code{diagEntry}. Note that most WGCNA functions expect the diagonal of an adjacency matrix to be 1. 
}
\value{
  The adjacency matrix that encodes the network.
}
\author{
Peter Langfelder
}

\seealso{
\code{adjacency} for calculation of a correlation network (adjacency) from a numeric matrix such as
expression data 

\code{adjacency.fromSimilarity} for simpler calculation of a network from a symmetric similarity matrix.
}
\keyword{misc}% __ONLY ONE__ keyword per line
