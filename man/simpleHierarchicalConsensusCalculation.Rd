\name{simpleHierarchicalConsensusCalculation}
\alias{simpleHierarchicalConsensusCalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple hierarchical consensus calculation
}
\description{
Hierarchical consensus calculation without calibration. 
}
\usage{
simpleHierarchicalConsensusCalculation(individualData, consensusTree, level = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individualData}{
Individual data from which the consensus is to be calculated. It can be either a list or a
\code{\link{multiData}} structure. Each element in \code{individulData} should be a numeric
object (vector, matrix or array).
}
  \item{consensusTree}{
A list specifying the consensus calculation. See details.
}
  \item{level}{
Integer which the user should leave at 1.  This serves to keep default set names unique.
}
}
\details{
This function calculates consensus in a hierarchical manner, using a separate (and possibly different) set of
consensus options at each step. The "recipe" for the consensus calculation is supplied in the argument
\code{consensusTree}.

The argument \code{consensusTree} should have the following components: (1) \code{inputs} must be either a
character vector whose components match \code{names(inputData)}, or consensus trees in the own right.
(2) \code{consensusOptions} must be a list of class \code{"ConsensusOptions"} that specifies options for
calculating the consensus. A suitable set of options can be obtained by calling
\code{\link{newConsensusOptions}}. (3) Optionally, the component \code{analysisName} can be a single
character string giving the name for the analysis. When intermediate results are returned, they are returned
in a list whose names will be set from \code{analysisName} components, if they exist.

Unlike the similar function \code{\link{hierarchicalConsensusCalculation}}, this function ignores the
calibration settings in the \code{consensusOptions} component of \code{consensusTree}; no calibration of
input data is performed.

The actual consensus calculation at each level of the consensus tree
is carried out in function \code{\link{simpleConsensusCalculation}}. The consensus options for each individual
consensus calculation are independent from one another, i.e., the consensus options for different steps can
be different.

}
\value{
A list with a single component  \code{consensus}, containing the consensus data of the same dimensions as the
individual entries in the input \code{individualData}. This perhaps somewhat cumbersome convention is used to
make the output compatible with that of \code{\link{hierarchicalConsensusCalculation}}.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{simpleConsensusCalculation}} for a "single-level" consensus calculation;

\code{\link{hierarchicalConsensusCalculation}} for hierarchical consensus calculation with calibration

}
\keyword{misc}% __ONLY ONE__ keyword per line
