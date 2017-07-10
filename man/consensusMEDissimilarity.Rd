\name{consensusMEDissimilarity}
\alias{consensusMEDissimilarity}
\title{ Consensus dissimilarity of module eigengenes. }
\description{
  Calculates consensus dissimilarity \code{(1-cor)} of given module eigengenes realized in several sets.
}
\usage{
consensusMEDissimilarity(MEs, useAbs = FALSE, useSets = NULL, method = "consensus")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MEs}{Module eigengenes of the same modules in several sets. }
  \item{useAbs}{Controls whether absolute value of correlation should be used instead of correlation in
the calculation of dissimilarity. }
  \item{useSets}{If the consensus is to include only a selection of the given sets, this vector (or
scalar in the case of a single set) can be used to specify the selection. If \code{NULL}, all sets will
be used. }
  \item{method}{A character string giving the method to use. Allowed values are (abbreviations of) 
\code{"consensus"} and \code{"majority"}. The consensus dissimilarity is calculated as the
minimum of given set dissimilarities for \code{"consensus"} and as the average for \code{"majority"}.}
}
\details{
This function calculates the
individual set dissimilarities of the given eigengenes in each set, then takes the (parallel) maximum or
average over all sets. For details on the structure of imput data, see \code{\link{checkSets}}. 
}
\value{
A dataframe containing the matrix of dissimilarities, with \code{names} and \code{rownames} set
appropriately.
}
%\references{
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{ \code{\link{checkSets}}}
%\examples{
\keyword{misc}% __ONLY ONE__ keyword per line
