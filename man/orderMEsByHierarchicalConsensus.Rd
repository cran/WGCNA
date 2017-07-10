\name{orderMEsByHierarchicalConsensus}
\alias{orderMEsByHierarchicalConsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Order module eigengenes by their hierarchical consensus similarity
}
\description{
This function calculates a hiearchical consensus similarity of the input eigengenes, clusters the eigengenes
according to the similarity and returns the input module eigengenes ordered by the order of resulting
dendrogram.
}
\usage{
orderMEsByHierarchicalConsensus(
    MEs, 
    networkOptions, 
    consensusTree, 
    greyName = "ME0", 
    calibrate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MEs}{
Module eigengenes, or more generally, vectors, to be ordered, in a \code{\link{multiData}} format: A vector of
lists, one per set. Each set must contain a component \code{data} that contains the module eigenegens or
general vectors, with
rows corresponding to samples and columns to genes or probes. }
  \item{networkOptions}{
A single list of class \code{\link{NetworkOptions}} giving options for network calculation for all of the
networks, or a \code{\link{multiData}} structure containing one such list for each input data set.
}
  \item{consensusTree}{
A list specifying the consensus calculation. See \code{\link{newConsensusTree}} for details.
}
  \item{greyName}{
Specifies the column name of eigengene of the "module" that contains unassigned genes. This eigengene
(column) will be excluded from the clustering and will be put last in the order.
}
  \item{calibrate}{
Logical: should module eigengene similarities be calibrated? This setting overrides the calibration options
in \code{consensusTree}.
}
}
\value{A \code{\link{multiData}} structure of the same format as the input \code{MEs}, with columns ordered
by the calculated dendrogram.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{hierarchicalConsensusMEDissimilarity}} for calculating the consensus ME dissimilarity
}
\keyword{misc}% __ONLY ONE__ keyword per line
