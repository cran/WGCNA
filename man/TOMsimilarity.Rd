\name{TOMsimilarity}
\alias{TOMsimilarity}
\alias{TOMdist}
\title{ Topological overlap matrix similarity and dissimilarity}
\description{
  Calculation of the topological overlap matrix, and the corresponding dissimilarity, 
from a given adjacency matrix.
}
\usage{
TOMsimilarity(
    adjMat,
    TOMType = "unsigned",
    TOMDenom = "min",
    suppressTOMForZeroAdjacencies = FALSE,
    useInternalMatrixAlgebra = FALSE,
    verbose = 1,
    indent = 0)
TOMdist(
    adjMat,
    TOMType = "unsigned",
    TOMDenom = "min",
    suppressTOMForZeroAdjacencies = FALSE,
    useInternalMatrixAlgebra = FALSE,
    verbose = 1,
    indent = 0)
}
\arguments{
  \item{adjMat}{ adjacency matrix, that is a square, symmetric matrix with entries between 0 and 1 
(negative values are allowed if \code{TOMType=="signed"}).  }

  \item{TOMType}{ a character string specifying TOM type to be calculated. One of \code{"unsigned"},
\code{"signed"}. If \code{"unsigned"}, the standard TOM will be used (more generally, TOM
function will receive the adjacency as input). If \code{"signed"}, TOM will keep track of the sign of
the adjacency between neighbors. }

  \item{TOMDenom}{ a character string specifying the TOM variant to be used. Recognized values are
\code{"min"} giving the standard TOM described in Zhang and Horvath (2005), and \code{"mean"} in which
the \code{min} function in the denominator is replaced by \code{mean}. The \code{"mean"} may produce 
better results but at this time should be considered experimental.}
%The default mean denominator variant
%is preferrable and we recommend using it unless the user needs to reproduce older results obtained using
%the standard, minimum denominator TOM. }

  \item{suppressTOMForZeroAdjacencies}{Logical: should TOM be set to zero for zero adjacencies?}

  \item{useInternalMatrixAlgebra}{Logical: should WGCNA's own, slow, matrix multiplication be used instead of
R-wide BLAS? Only useful for debugging.}

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{

  The functions perform basically the same calculations of topological overlap. \code{TOMdist} turns the
overlap (which is a measure of similarity) into a measure of dissimilarity by subtracting it from 1. 

  Basic checks on the adjacency matrix are performed and missing entries are replaced by zeros.
  If \code{TOMType = "unsigned"}, entries of the adjacency matrix are required to lie between 0 and 1;
for \code{TOMType = "signed"} they can be between -1 and 1. In both cases the resulting TOM entries, as
well as the corresponding dissimilarities, lie between 0 and 1. 

  The underlying C code assumes that the diagonal of the adjacency matrix equals 1. If this is not the case,
the diagonal of the input is set to 1 before the calculation begins.

}
\value{
  A matrix holding the topological overlap.
}
\references{ Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
}
\author{ Peter Langfelder }
\seealso{ \code{\link{TOMsimilarityFromExpr}} }
\keyword{ misc }% __ONLY ONE__ keyword per line
