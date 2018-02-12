\name{subsetTOM}
\alias{subsetTOM}
\title{ Topological overlap for a subset of a whole set of genes }
\description{
  This function calculates topological overlap of a subset of vectors with respect to a whole data
set.
}
\usage{
subsetTOM(
  datExpr, 
  subset,
  corFnc = "cor", corOptions = "use = 'p'", 
  weights = NULL,  
  networkType = "unsigned", 
  power = 6, 
  verbose = 1, indent = 0)
}
\arguments{

  \item{datExpr}{ a data frame containing the expression data of the whole set, 
with rows corresponding to samples and columns to genes. }

  \item{subset}{ a single logical or numeric vector giving the indices of the nodes for which the TOM is to
be calculated. }

  \item{corFnc}{ character string giving the correlation function to be used for the adjacency
calculation. Recommended choices are \code{"cor"} and \code{"bicor"}, but other functions can be used as
well. }

  \item{corOptions}{  character string giving further options to be passed to the correlation function. }

  \item{weights}{optional observation weights for \code{datExpr} to be used in correlation calculation.
A matrix of the same dimensions as \code{datExpr}, containing non-negative weights. Only used with Pearson
correlation.}


  \item{networkType}{ character string giving network type. Allowed values are (unique abbreviations of)
\code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{power}{ soft-thresholding power for network construction. }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
This function is designed to calculated topological overlaps of small subsets of large expression data sets,
for example in individual modules.
}
\value{
  A matrix of dimensions \code{n*n}, where \code{n} is the number of entries selected by \code{block}.
}
\references{ Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
}

\author{ Peter Langfelder }
\seealso{ \code{\link{TOMsimilarity}} for standard calculation of topological overlap. }
\keyword{ misc }
