\name{hierarchicalConsensusTOM}
\alias{hierarchicalConsensusTOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of hierarchical consensus topological overlap matrix
}
\description{
This function calculates consensus topological overlap in a hierarchical manner.
}
\usage{
hierarchicalConsensusTOM(
      # ... information needed to calculate individual TOMs
      multiExpr,
      multiWeights = NULL,

      # Data checking options
      checkMissingData = TRUE,

      # Blocking options
      blocks = NULL,
      maxBlockSize = 20000,
      blockSizePenaltyPower = 5,
      nPreclusteringCenters = NULL,
      randomSeed = 12345,

      # Network construction options
      networkOptions,

      # Save individual TOMs?

      keepIndividualTOMs = TRUE,
      individualTOMFileNames = "individualTOM-Set\%s-Block\%b.RData",

      # ... or information about individual (more precisely, input) TOMs
      individualTOMInfo = NULL,

      # Consensus calculation options 
      consensusTree,

      useBlocks = NULL,

      # Save calibrated TOMs?
      saveCalibratedIndividualTOMs = FALSE,
      calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set\%s-Block\%b.RData",

      # Return options
      saveConsensusTOM = TRUE,
      consensusTOMFilePattern = "consensusTOM-\%a-Block\%b.RData",
      getCalibrationSamples = FALSE,

      # Return the intermediate results as well?  
      keepIntermediateResults = saveConsensusTOM,

      # Internal handling of TOMs
      useDiskCache = NULL, 
      chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      # Behavior
      collectGarbage = TRUE,
      verbose = 1,
      indent = 0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{ Expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{multiWeights}{ optional observation weights in the same format (and dimensions) as \code{multiExpr}.
These weights are used for correlation calculations with data in  \code{multiExpr}.}

  \item{checkMissingData}{
Logical: should data be checked for excessive numbers of missing entries in
genes and samples, and for genes with zero variance? See details.
}
  \item{blocks}{
Optional specification of blocks in which hierarchical clustering and module detection
should be performed. If given, must be a numeric vector with one entry per gene
of \code{multiExpr} giving the number of the block to which the corresponding gene belongs.
}
  \item{maxBlockSize}{Integer giving maximum block size for module detection. Ignored if \code{blocks}
above is non-NULL. Otherwise, if the number of genes in \code{datExpr} exceeds \code{maxBlockSize}, genes
will be pre-clustered into blocks whose size should not exceed \code{maxBlockSize}. }

   \item{blockSizePenaltyPower}{Number specifying how strongly blocks should be penalized for exceeding the
maximum size. Set to a lrge number or \code{Inf} if not exceeding maximum block size is very important.}

 \item{nPreclusteringCenters}{Number of centers to be used in the preclustering. Defaults to smaller of
\code{nGenes/20} and \code{100*nGenes/maxBlockSize}, where \code{nGenes} is the nunber of genes (variables)
in \code{multiExpr}.}

 \item{randomSeed}{Integer to be used as seed for the random number generator before the function
starts. If a current seed exists, it is saved and restored upon exit. If \code{NULL} is given, the
function will not save and restore the seed. }

  \item{networkOptions}{
A single list of class \code{\link{NetworkOptions}} giving options for network calculation for all of the
networks, or a \code{\link{multiData}} structure containing one such list for each input data set.
}
  \item{keepIndividualTOMs}{
Logical: should individual TOMs be retained after the calculation is finished?
}
  \item{individualTOMFileNames}{
Character string giving the file names to save individual TOMs into. The
following tags should be used to make the file names unique for each set and block: \code{\%s} will be
replaced by the set number; \code{\%N} will be replaced by the set name (taken from \code{names(multiExpr)})
if it exists, otherwise by set number; \code{\%b} will be replaced by the block number. If the file names
turn out to be non-unique, an error will be generated.}

  \item{individualTOMInfo}{
A list, typically returned by \code{\link{individualTOMs}}, containing information about the topological
overlap matrices in the individual data sets in \code{multiExpr}. See the output of 
\code{\link{individualTOMs}} for
details on the content of the list.
}

  \item{consensusTree}{
A list specifying the consensus calculation. See details.
}
  \item{useBlocks}{
Optional vector giving the blocks that should be used for the calcualtions. If \code{NULL}, all 
all blocks will be used.
}
  \item{saveCalibratedIndividualTOMs}{
Logical: should the calibrated individual TOMs be saved?
}
  \item{calibratedIndividualTOMFilePattern}{
Specification of file names in which calibrated individual TOMs should be saved.
}
  \item{saveConsensusTOM}{
Logical: should the consensus TOM be saved to disk?
}
  \item{consensusTOMFilePattern}{
Character string giving the file names to save consensus TOMs into. The
following tags should be used to make the file names unique for each set and block: \code{\%s} will be
replaced by the set number; \code{\%N} will be replaced by the set name (taken from \code{names(multiExpr)})
if it exists, otherwise by set number; \code{\%b} will be replaced by the block number. If the file names
turn out to be non-unique, an error will be generated.
}
  \item{getCalibrationSamples}{
  Logical: should the sampled values used for network calibration be returned?
}
  \item{keepIntermediateResults}{
Logical: should intermediate consensus TOMs be saved as well?
}
  \item{useDiskCache}{
Logical: should disk cache be used for consensus calculations? The disk cache can be used to store chunks of
calibrated data that are small enough to fit one chunk from each set into memory (blocks may be small enough
to fit one block of one set into memory, but not small enough to fit one block from all sets in a consensus
calculation into memory at the same time). Using disk cache is slower but lessens the memory footprint of
the calculation.
As a general guide, if individual data are split into blocks, we
recommend setting this argument to \code{TRUE}. If this argument is \code{NULL}, the function will decide
whether to use disk cache based on the number of sets and block sizes.
}
  \item{chunkSize}{ network similarities are saved in smaller chunks of size \code{chunkSize}. If \code{NULL},
an appropriate chunk size will be determined from an estimate of available memory. Note that if the chunk size
is greater than the memory required for storing intemediate results, disk cache use will automatically be
disabled. }

  \item{cacheDir}{ character string containing the directory into which cache files should be written. The
user should make sure that the filesystem has enough free space to hold the cache files which can get quite
large.
}
  \item{cacheBase}{ character string containing the desired name for the cache files. The actual file
names will consists of \code{cacheBase} and a suffix to make the file names unique. }

  \item{collectGarbage}{
Logical: should garbage be collected after memory-intensive operations?
}
  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
This function is essentially a wrapper for \code{\link{hierarchicalConsensusCalculation}}, with a few
additional operations specific to calculations of topological overlaps.
}
\value{
A list that contains the output of \code{\link{hierarchicalConsensusCalculation}} and two extra components:
\item{individualTOMInfo}{A copy of the input \code{individualTOMInfo} if it was non-\code{NULL}, or the
result of \code{\link{individualTOMs}}.
}
\item{consensusTree}{A copy of the input \code{consensusTree}.}
}

\author{
Peter Langfelder
}
\seealso{
\code{\link{hierarchicalConsensusCalculation}} for the actual hierarchical consensus calculation;

\code{\link{individualTOMs}} for the calculation of individual TOMs in a format suitable for consensus
calculation.
}
\keyword{misc }% __ONLY ONE__ keyword per line
