\name{unsignedAdjacency}
\alias{unsignedAdjacency}
\title{ Calculation of unsigned adjacency }
\description{
  Calculation of the unsigned network adjacency from expression data. The restricted set of parameters
for this function should allow a faster and less memory-hungry calculation.
}
\usage{
unsignedAdjacency(
  datExpr, 
  datExpr2 = NULL, 
  power = 6, 
  corFnc = "cor", corOptions = "use = 'p'")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{ expression data. A data frame in which columns are genes and rows ar
samples. Missing values are ignored. }

  \item{datExpr2}{ optional specification of a second set of expression data. See details. }

  \item{power}{ soft-thresholding power for network construction. }

  \item{corFnc}{ character string giving the correlation function to be used for the adjacency
calculation. Recommended choices are \code{"cor"} and \code{"bicor"}, but other functions can be used as
well. }

  \item{corOptions}{  character string giving further options to be passed to the correlation function }

}
\details{
  The correlation function will be called with arguments \code{datExpr, datExpr2} plus any extra
arguments given in \code{corOptions}. If \code{datExpr2} is \code{NULL}, 
the standard correlation functions will calculate the corelation of columns in \code{datExpr}.
}
\value{
  Adjacency matrix of dimensions \code{n*n}, where \code{n} is the number of genes in \code{datExpr}.
}
\references{
  Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
 }
\author{ Steve Horvath and Peter Langfelder }
\seealso{ \code{\link{adjacency}} }
\keyword{ misc }% __ONLY ONE__ keyword per line
