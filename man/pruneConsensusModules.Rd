\name{pruneConsensusModules}
\alias{pruneConsensusModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prune (hierarchical) consensus modules by removing genes with low eigengene-based intramodular connectivity
}
\description{
This function prunes (hierarchical) consensus modules by removing genes with low eigengene-based intramodular
connectivity (KME) and by removing modules that do not have a certain minimum number of genes with a required
minimum KME.
}
\usage{
pruneConsensusModules(  multiExpr,
  multiWeights = NULL,
  multiExpr.imputed = NULL,
  MEs = NULL,
  labels,

  unassignedLabel = if (is.numeric(labels)) 0 else "grey",

  networkOptions,
  consensusTree,

  minModuleSize,
  minCoreKMESize = minModuleSize/3,
  minCoreKME = 0.5, 
  minKMEtoStay = 0.2,

  # Module eigengene calculation options
  impute = TRUE, 
  collectGarbage = FALSE,
  checkWeights = TRUE,

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
\item{MEs}{Optional consensus module eigengenes, in multi-set format analogous to that of \code{multiExpr}.}
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

  \item{collectGarbage}{
Logical: should garbage be collected after some of the memory-intensive steps?
}

  \item{checkWeights}{Logical: should \code{multiWeights} be checked to make sure their dimensions are
concordant with \code{multiExpr} and the weights are valid?}

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\value{
The pruned module labels: a vector of the same form as the input \code{labels}.
}
\author{
Peter Langfelder
}
\keyword{misc}
