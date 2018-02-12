\name{pruneAndMergeConsensusModules}
\alias{pruneAndMergeConsensusModules}
\title{
Iterative pruning and merging of (hierarchical) consensus modules
}
\description{
This function prunes genes with low consensus eigengene-based intramodular connectivity (kME) from modules
and merges modules whose consensus similarity is high. The process is repeated until the modules become
stable.
}
\usage{
pruneAndMergeConsensusModules(
  multiExpr,
  multiWeights = NULL,
  multiExpr.imputed = NULL,
  labels,
  unassignedLabel = if (is.numeric(labels)) 0 else "grey",
  networkOptions,
  consensusTree,

  # Pruning options
  minModuleSize,
  minCoreKMESize = minModuleSize/3,
  minCoreKME = 0.5, 
  minKMEtoStay = 0.2,

  # Module eigengene calculation and merging options
  impute = TRUE,
  trapErrors = FALSE,
  calibrateMergingSimilarities = FALSE,
  mergeCutHeight = 0.15,

  # Behavior
  iterate = TRUE,
  collectGarbage = FALSE,
  getDetails = TRUE,
  verbose = 1, indent=0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{ Expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{multiWeights}{ optional observation weights in the same format (and dimensions) as \code{multiExpr}.
These weights are used for correlation calculations with data in  \code{multiExpr}.}

  \item{multiExpr.imputed}{If \code{multiExpr} contain missing data, this argument can be used to supply the
expression data with missing data imputed. If not given, the \code{\link[impute]{impute.knn}} function will
be used to impute the missing data.}

\item{labels}{
A vector (numeric, character or a factor) giving module labels for each variable (gene) in multiExpr.
}
  \item{unassignedLabel}{
The label (value in \code{labels})
that represents unassigned genes. Module of this label will
not enter the module eigengene clustering and will not be merged with other modules.}
  \item{networkOptions}{
A single list of class \code{\link{NetworkOptions}} giving options for network calculation for all of the
networks, or a \code{\link{multiData}} structure containing one such list for each input data set.
}
  \item{consensusTree}{
A list of class \code{\link{ConsensusTree}} specifying the consensus calculation.
}

  \item{minModuleSize}{Minimum number of genes in a module. Modules that have fewer genes (after trimming)
will be removed (i.e., their genes will be given the unassigned label).}

  \item{minCoreKME}{ a number between 0 and 1. If a detected module does not have at least
\code{minModuleKMESize} genes with consensus eigengene connectivity at least \code{minCoreKME}, the module is
disbanded (its genes are unlabeled).}

  \item{minCoreKMESize}{ see \code{minCoreKME} above. }

  \item{minKMEtoStay}{ genes whose consensus eigengene connectivity to their module eigengene is lower than
\code{minKMEtoStay} are removed from the module.}

  \item{impute}{ logical: should imputation be used for module eigengene calculation? See
\code{\link{moduleEigengenes}} for more details. }

 \item{trapErrors}{ logical: should errors in calculations be trapped? }
  \item{calibrateMergingSimilarities}{
Logical: should module eigengene similarities be calibrated before calculating the consensus? Although
calibration is in principle desirable, the calibration methods currently available assume large data and do
not work very well on eigengene similarities.
}
  \item{mergeCutHeight}{
Dendrogram cut height for module merging.
}
  \item{iterate}{
Logical: should the pruning and merging process be iterated until no changes occur? If \code{FALSE}, only one
iteration will be carried out.
}
  \item{collectGarbage}{
Logical: should garbage be collected after some of the memory-intensive steps?
}
  \item{getDetails}{
Logical: should certain intermediate results be returned? These include labels and module merging information
at each iteration (see return value).
}
  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
%\details{
%Each iteration 
%}
\value{
If input \code{getDetails} is \code{FALSE}, a vector the resulting module labels. If \code{getDetails} is
\code{TRUE}, a list with these components:
\item{labels}{The resulting module labels}
\item{details}{A list. The first component, named \code{originalLabels}, contains a copy of the input labels.
The following components are named \code{Iteration.1}, \code{Iteration.2} etc and contain, for each iteration,
components \code{prunedLabels} (the result of pruning in that iteration) and \code{mergeInfo} (result of the
call to \code{\link{hierarchicalMergeCloseModules}} in that iteration).}
}
\author{
Peter Langfelder
}
\seealso{
The underlying functions \code{\link{pruneConsensusModules}} and \code{\link{hierarchicalMergeCloseModules}}.
}
\keyword{misc}
