\name{individualTOMs}
\alias{individualTOMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate individual correlation network matrices
}
\description{
This function calculates correlation network matrices (adjacencies or topological overlaps), after optionally
first pre-clustering input data into blocks.
}
\usage{
individualTOMs(
   multiExpr,
   multiWeights = NULL,
   multiExpr.imputed = NULL,  

   # Data checking options
   checkMissingData = TRUE,

   # Blocking options
   blocks = NULL,
   maxBlockSize = 5000,
   blockSizePenaltyPower = 5,
   nPreclusteringCenters = NULL,
   randomSeed = 54321,

   # Network construction options
   networkOptions,

   # Save individual TOMs? 
   saveTOMs = TRUE,
   individualTOMFileNames = "individualTOM-Set\%s-Block\%b.RData",

   # Behaviour options
   collectGarbage = TRUE,
   verbose = 2, indent = 0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{ expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{multiWeights}{ optional observation weights in the same format (and dimensions) as \code{multiExpr}.
These weights are used for correlation calculations with data in  \code{multiExpr}.}

  \item{multiExpr.imputed}{
Optional version of \code{multiExpr} with missing data imputed. If not given and \code{multiExpr} contains
missing data, they will be imputed using the function \code{\link[impute]{impute.knn}}.
}
  \item{checkMissingData}{logical: should data be checked for excessive numbers of missing entries in
genes and samples, and for genes with zero variance? See details. }

    \item{blocks}{ optional specification of blocks in which hierarchical clustering and module detection
should be performed. If given, must be a numeric vector with one entry per gene
of \code{multiExpr} giving the number of the block to which the corresponding gene belongs. }

  \item{maxBlockSize}{ integer giving maximum block size for module detection. Ignored if \code{blocks}
above is non-NULL. Otherwise, if the number of genes in \code{datExpr} exceeds \code{maxBlockSize}, genes
will be pre-clustered into blocks whose size should not exceed \code{maxBlockSize}. }

   \item{blockSizePenaltyPower}{number specifying how strongly blocks should be penalized for exceeding the
maximum size. Set to a lrge number or \code{Inf} if not exceeding maximum block size is very important.}

 \item{nPreclusteringCenters}{number of centers to be used in the preclustering. Defaults to smaller of
\code{nGenes/20} and \code{100*nGenes/maxBlockSize}, where \code{nGenes} is the nunber of genes (variables)
in
\code{multiExpr}.}

 \item{randomSeed}{ integer to be used as seed for the random number generator before the function
starts. If a current seed exists, it is saved and restored upon exit. If \code{NULL} is given, the
function will not save and restore the seed. }

  \item{networkOptions}{
A single list of class \code{\link{NetworkOptions}} giving options for network calculation for all of the
networks, or a \code{\link{multiData}} structure containing one such list for each input data set. 
}
  \item{saveTOMs}{logical: should individual TOMs be saved to disk (\code{TRUE}) or retuned directly in the
return value (\code{FALSE})?}

  \item{individualTOMFileNames}{character string giving the file names to save individual TOMs into. The
following tags should be used to make the file names unique for each set and block: \code{\%s} will be
replaced by the set number; \code{\%N} will be replaced by the set name (taken from \code{names(multiExpr)})
if it exists, otherwise by set number; \code{\%b} will be replaced by the block number. If the file names
turn out to be non-unique, an error will be generated.}

  \item{collectGarbage}{
Logical: should garbage collection be called after each block calculation? This can be useful when the data
are large, but could unnecessarily slow down calculation with small data.
}
  \item{verbose}{ Integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ Indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
 The function starts by optionally filtering out samples that have too many missing entries and genes
that have either too many missing entries or zero variance in at least one set. Genes that are filtered
out are excluded from the network calculations.

If \code{blocks} is not given and
the number of genes (columns) in \code{multiExpr} 
exceeds \code{maxBlockSize}, genes are pre-clustered into blocks using the function
\code{\link{consensusProjectiveKMeans}}; otherwise all genes are treated in a single block. Any missing data
in \code{multiExpr} will be imputed; if imputed data are already available, they can be supplied separately.

For each block of genes, the network adjacency is constructed and (if requested) topological overlap is calculated
in each set. The topological overlaps can be saved to disk as RData files, or returned directly within the
return value (see below). Note that the matrices can be big and returning them within the return value can
quickly exhaust the system's memory. In particular, if the block-wise calculation is necessary, it is
usually impossible to return all matrices in the return value.
}
\value{
A list with the following components:
\item{blockwiseAdjacencies}{A \code{\link{multiData}} structure containing (possibly blockwise) network
matrices for each input data set. The network matrices are stored as \code{\link{BlockwiseData}} objects.}
\item{setNames}{A copy of \code{names(multiExpr)}.}
\item{nSets}{Number of sets in \code{multiExpr}}
\item{blockInfo}{A list of class \code{\link{BlockInformation}}, giving information about blocks and gene and
sample filtering.}
\item{networkOptions}{The input \code{networkOptions}, returned as a \code{\link{multiData}} structure with
one entry per input data set.}
}
\author{
Peter Langfelder
}
\seealso{
Input arguments and output components of this function use \code{\link{multiData}},
\code{\link{NetworkOptions}}, \code{\link{BlockwiseData}}, and \code{\link{BlockInformation}}.

Underlying functions of interest include \code{\link{consensusProjectiveKMeans}},
\code{\link{TOMsimilarityFromExpr}}.
}
\keyword{misc}
