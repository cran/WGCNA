\name{newBlockInformation}
\alias{newBlockInformation}
\alias{BlockInformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a list holding information about dividing data into blocks
}
\description{
This function creates a list storing information about dividing data into blocks, as well as about possibly
excluding genes or samples with excessive numbers of missing data.
}
\usage{
newBlockInformation(blocks, goodSamplesAndGenes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blocks}{
A vector giving block labels. It is assumed to be a numeric vector with block labels consecutive integers
starting at 1. 
}
  \item{goodSamplesAndGenes}{
A list returned by \code{\link{goodSamplesGenes}} or \code{\link{goodSamplesGenesMS}}.
}
}
\value{
A list with \code{class} attribute set to \code{BlockInformation}, with the following componens:
\item{blocks}{A copy of the input \code{blocks}.}
\item{blockGenes}{A list with one component per block, giving the indices of elements in \code{block} whose
value is the same.}
\item{goodSamplesAndGenes}{A copy of input \code{goodSamplesAndGenes}.}
\item{nGGenes}{Number of `good' genes in \code{goodSamplesAndGenes}.}
\item{gBlocks}{The input \code{blocks} restricted to `good' genes in \code{goodSamplesAndGenes}.}
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{goodSamplesGenes}}, \code{\link{goodSamplesGenesMS}}.
}
\keyword{misc}
