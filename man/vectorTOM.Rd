\name{vectorTOM}
\alias{vectorTOM}
\title{ Topological overlap for a subset of the whole set of genes }
\description{
  This function calculates topological overlap of a small set of vectors with respect to a whole data
set.
}
\usage{
vectorTOM(
  datExpr, 
  vect, 
  subtract1 = FALSE, 
  blockSize = 2000, 
  corFnc = "cor", corOptions = "use = 'p'", 
  networkType = "unsigned", 
  power = 6, 
  verbose = 1, indent = 0)
}
\arguments{

  \item{datExpr}{ a data frame containing the expression data of the whole set, 
with rows corresponding to samples and columns to genes. }

  \item{vect}{ a single vector or a matrix-like object containing vectors whose topological overlap is to
be calculated. }

  \item{subtract1}{ logical: should calculation be corrected for self-correlation? Set this to
\code{TRUE} if \code{vect} contains a subset of \code{datExpr}. }

  \item{blockSize}{ maximum block size for correlation calculations. Only important if \code{vect}
contains a large number of columns. }

  \item{corFnc}{ character string giving the correlation function to be used for the adjacency
calculation. Recommended choices are \code{"cor"} and \code{"bicor"}, but other functions can be used as
well. }

  \item{corOptions}{  character string giving further options to be passed to the correlation function. }

  \item{networkType}{ character string giving network type. Allowed values are (unique abbreviations of)
\code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{power}{ soft-thresholding power for network construction. }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
  Topological overlap can be viewed as the normalized count of shared neighbors encoded in an adjacency
matrix. In this case, the adjacency matrix is calculated between the columns of \code{vect} and
\code{datExpr} and the topological overlap of vectors in \code{vect} measures the number of shared
neighbors in \code{datExpr} that vectors of \code{vect} share. 
}
\value{
  A matrix of dimensions \code{n*n}, where \code{n} is the number of columns in \code{vect}.
}
\references{ Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
}

\author{ Peter Langfelder }
\seealso{ \code{\link{TOMsimilarity}} for standard calculation of topological overlap. }
\keyword{ misc }
