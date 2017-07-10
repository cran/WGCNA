\name{hierarchicalConsensusMEDissimilarity}
\alias{hierarchicalConsensusMEDissimilarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical consensus calculation of module eigengene dissimilarity
}
\description{
Hierarchical consensus calculation of module eigengene dissimilarities, or more generally, correlation-based
dissimilarities of sets of vectors.
}
\usage{
hierarchicalConsensusMEDissimilarity(
   MEs, 
   networkOptions, 
   consensusTree, 
   greyName = "ME0", 
   calibrate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MEs}{
A \code{\link{multiData}} structure containing vectors (usually module eigengenes) whose consensus
dissimilarity is to be calculated.
}
  \item{networkOptions}{
A \code{\link{multiData}} structure containing, for each input data set, a list of class \code{\link{NetworkOptions}} giving options for network calculation for all of the networks.
}
  \item{consensusTree}{
A list specifying the consensus calculation. See details.
}
  \item{greyName}{
Name of the "grey" module eigengene. Currently not used.
}
  \item{calibrate}{
Logical: should the dissimilarities be calibrated using the calibration method specified in
\code{consensusTree}? See details.
}
}
\details{
This function first calculates the similarities of the ME vectors from their correlations, using the appropriate
options in \code{networkOptions} (correlation type and options, signed or unsigned dissimilarity etc). This
results in a similarity matrix in each of the input data sets. 

Next, a hierarchical consensus of the similarities is calculated via a call to
\code{\link{hierarchicalConsensusCalculation}}, using the consensus specification and
options in \code{consensusTree}. In typical use, \code{consensusTree} contains the same consensus
specification as the consensus network calculation that gave rise to the consensus modules whose eigengenes
are contained in \code{MEs} but this is not mandatory.

The argument \code{consensusTree} should have the following components: (1) \code{inputs} must be either a
character vector whose components match \code{names(inputData)}, or consensus trees in the own right.
(2) \code{consensusOptions} must be a list of class \code{"ConsensusOptions"} that specifies options for
calculating the consensus. A suitable set of options can be obtained by calling
\code{\link{newConsensusOptions}}. (3) Optionally, the component \code{analysisName} can be a single
character string giving the name for the analysis. When intermediate results are returned, they are returned
in a list whose names will be set from \code{analysisName} components, if they exist.

In the final step, the consensus similarity is turned into a dissimilarity by subtracting it from 1.
}
\value{
A matrix with rows and columns corresponding to the variables (modules) in MEs, containing the consensus
dissimilarities. 
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{hierarchicalConsensusCalculation}} for the actual consensus calculation.
}
\keyword{misc}
