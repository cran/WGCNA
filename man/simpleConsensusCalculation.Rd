\name{simpleConsensusCalculation}
\alias{simpleConsensusCalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple calculation of a single consenus
}
\description{
This function calculates a single consensus from given individual data.
}
\usage{
simpleConsensusCalculation(
  individualData, 
  consensusOptions, 
  verbose = 1, 
  indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individualData}{
Individual data from which the consensus is to be calculated. It can be either a list or a
\code{\link{multiData}} structure in which each element is a numeric vector or array.
}
  \item{consensusOptions}{
A list of class \code{ConsensusOptions} that contains options for the consensus calculation. A suitable list
can be obtained by calling function \code{\link{newConsensusOptions}}.
}
  \item{verbose}{Integer level of verbosity of diagnostic messages. Zero means silent, higher values make the
output
progressively more and more verbose. }

  \item{indent}{Indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
Consensus is defined as the element-wise (also known as "parallel") quantile of of the individual data at
probability given by the \code{consensusQuantile} element of \code{consensusOptions}. 
}
\value{
A numeric vector or array of the same dimensions as each element of \code{individualData}
}
\references{
Consensus network analysis was originally described in
Langfelder P, Horvath S. Eigengene networks for studying the relationships
between co-expression modules. BMC Systems Biology 2007, 1:54
http://www.biomedcentral.com/1752-0509/1/54
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{consensusCalculation}} for consensus calculation that can work with \code{\link{BlockwiseData}}
and can calibrate data before calculating consensus.
}
\keyword{misc}
