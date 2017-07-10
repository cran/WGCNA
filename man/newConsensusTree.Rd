\name{newConsensusTree}
\alias{newConsensusTree}
\alias{ConsensusTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a new consensus tree
}
\description{
This function creates a new consensus tree, a class for representing "recipes" for hierarchical consensus
calculations.
}
\usage{
newConsensusTree(
   consensusOptions = newConsensusOptions(), 
   inputs, 
   analysisName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{consensusOptions}{
An object of class \code{ConsensusOptions}, usually obtained by calling \code{\link{newConsensusOptions}}.
}
  \item{inputs}{
A vector (or list) of inputs. Each component can be either a character string giving a names of a data set,
or another \code{ConsensusTree} object.
}
  \item{analysisName}{
Optional specification of a name for this consensus analysis. While this has no effect on the actual
consensus calculation, some functions use this character string to make certain file names unique.
}
}
\details{
Consensus trees specify a "recipe" for the calculation of hierarchical consensus in
\code{\link{hierarchicalConsensusCalculation}} and other functions. 
}
\value{
A list with class set to \code{"ConsensusTree"} with these components:
\code{consensusOptions}{A copy of the input \code{consensusOptions}.}
\code{inputs}{A copy of the input \code{inputs}.}
\code{analysisName}{A copy of the input \code{analysisName}.}
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{hierarchicalConsensusCalculation}} for hierarchical consensus calculation for which a
\code{ConsensusTree} object specifies the recipe
}
\keyword{misc}% __ONLY ONE__ keyword per line
