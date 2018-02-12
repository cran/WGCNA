\name{hierarchicalConsensusModules}
\alias{hierarchicalConsensusModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical consensus network construction and module identification
}
\description{
Hierarchical consensus network construction and module identification across multiple data sets.
}
\usage{
hierarchicalConsensusModules(
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
   randomSeed = 12345,

   # Network construction options. 
   networkOptions,

   # Save individual TOMs?
   saveIndividualTOMs = TRUE,
   individualTOMFileNames = "individualTOM-Set\%s-Block\%b.RData",
   keepIndividualTOMs = FALSE,

   # Consensus calculation options
   consensusTree = NULL,  

   # Return options
   saveConsensusTOM = TRUE,
   consensusTOMFilePattern = "consensusTOM-\%a-Block\%b.RData",

   # Keep the consensus? 
   keepConsensusTOM = saveConsensusTOM,

   # Internal handling of TOMs
   useDiskCache = NULL, chunkSize = NULL,
   cacheBase = ".blockConsModsCache",
   cacheDir = ".",

   # Alternative consensus TOM input from a previous calculation 
   consensusTOMInfo = NULL,

   # Basic tree cut options 
   deepSplit = 2, 
   detectCutHeight = 0.995, minModuleSize = 20,
   checkMinModuleSize = TRUE,

   # Advanced tree cut opyions
   maxCoreScatter = NULL, minGap = NULL,
   maxAbsCoreScatter = NULL, minAbsGap = NULL,
   minSplitHeight = NULL, minAbsSplitHeight = NULL,

   useBranchEigennodeDissim = FALSE,
   minBranchEigennodeDissim = mergeCutHeight,

   stabilityLabels = NULL,
   stabilityCriterion = c("Individual fraction", "Common fraction"),
   minStabilityDissim = NULL,

   pamStage = TRUE,  pamRespectsDendro = TRUE,

   iteratePruningAndMerging = FALSE,
   minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
   minKMEtoStay = 0.2,

   # Module eigengene calculation options

   impute = TRUE,
   trapErrors = FALSE,

   # Module merging options

   calibrateMergingSimilarities = FALSE,
   mergeCutHeight = 0.15, 
                    
   # General options
   collectGarbage = TRUE,
   verbose = 2, indent = 0,
   ...)
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
  \item{saveIndividualTOMs}{
Logical: should individual TOMs be saved to disk (\code{TRUE}) or retuned directly in the
return value (\code{FALSE})?
}
  \item{individualTOMFileNames}{
Character string giving the file names to save individual TOMs into. The
following tags should be used to make the file names unique for each set and block: \code{\%s} will be
replaced by the set number; \code{\%N} will be replaced by the set name (taken from \code{names(multiExpr)})
if it exists, otherwise by set number; \code{\%b} will be replaced by the block number. If the file names
turn out to be non-unique, an error will be generated.}
  \item{keepIndividualTOMs}{
Logical: should individual TOMs be retained after the calculation is finished?
}
  \item{consensusTree}{
A list specifying the consensus calculation. See details.
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
  \item{keepConsensusTOM}{
Logical: should consensus TOM be retained after the calculation ends? Depending on \code{saveConsensusTOM},
the retained TOM is either saved to disk or returned within the return value.
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
  \item{chunkSize}{
Integer giving the chunk size. If left \code{NULL}, a suitable size will be chosen automatically.
}
  \item{cacheDir}{
Directory in which to save cache files. The files are deleted on normal exit but persist if the function
terminates abnormally.
}
  \item{cacheBase}{
Base for the file names of cache files.
}
  \item{consensusTOMInfo}{
If the consensus TOM has been pre-calculated using function \code{\link{hierarchicalConsensusTOM}},
this argument can be used to supply it. If given, the consensus
TOM calculation options above are ignored. 
}
  \item{deepSplit}{Numeric value between 0 and 4. Provides a simplified control over how sensitive
module detection should be to module splitting, with 0 least and 4 most sensitive. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. 
}
  \item{detectCutHeight}{Dendrogram cut height for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  
}
  \item{minModuleSize}{Minimum module size for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  
}

  \item{checkMinModuleSize}{ logical: should sanity checks be performed on \code{minModuleSize}?}

  \item{maxCoreScatter}{ maximum scatter of the core for a branch to be a cluster, given as the fraction
of \code{cutHeight} relative to the 5th percentile of joining heights. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  }

  \item{minGap}{ minimum cluster gap given as the fraction of the difference between \code{cutHeight} and
the 5th percentile of joining heights. See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{maxAbsCoreScatter}{ maximum scatter of the core for a branch to be a cluster given as absolute
heights. If given, overrides \code{maxCoreScatter}. See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{minAbsGap}{ minimum cluster gap given as absolute height difference. If given, overrides
\code{minGap}. See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{minSplitHeight}{Minimum split height given as the fraction of the difference between
\code{cutHeight} and the 5th percentile of joining heights. Branches merging below this height will
automatically be merged. Defaults to zero but is used only if \code{minAbsSplitHeight} below is
\code{NULL}.}

  \item{minAbsSplitHeight}{Minimum split height given as an absolute height.
Branches merging below this height will automatically be merged. If not given (default), will be determined
from \code{minSplitHeight} above.}

   \item{useBranchEigennodeDissim}{Logical: should branch eigennode (eigengene) dissimilarity be considered
when merging branches in Dynamic Tree Cut?}

   \item{minBranchEigennodeDissim}{Minimum consensus branch eigennode (eigengene) dissimilarity for
branches to be considerd separate. The branch eigennode dissimilarity in individual sets
is simly 1-correlation of the
eigennodes; the consensus is defined as quantile with probability \code{consensusQuantile}.}

   \item{stabilityLabels}{Optional matrix of cluster labels that are to be used for calculating branch
dissimilarity based on split stability. The number of rows must equal the number of genes in
\code{multiExpr}; the number of columns (clusterings) is arbitrary. See
\code{\link{branchSplitFromStabilityLabels}} for details.}

  \item{stabilityCriterion}{One of \code{c("Individual fraction", "Common fraction")}, indicating which method
for assessing stability similarity of two branches should be used. We recommend \code{"Individual fraction"}
which appears to perform better; the \code{"Common fraction"} method is provided for backward compatibility
since it was the (only) method available prior to WGCNA version 1.60.}

   \item{minStabilityDissim}{Minimum stability dissimilarity criterion for two branches to be considered
separate. Should be a number between 0 (essentially no dissimilarity required) and 1 (perfect dissimilarity
or distinguishability based on \code{stabilityLabels}). See
\code{\link{branchSplitFromStabilityLabels}} for details.}

  \item{pamStage}{ logical.  If TRUE, the second (PAM-like) stage of module detection will be performed.
     See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{pamRespectsDendro}{Logical, only used when \code{pamStage} is \code{TRUE}.
If \code{TRUE}, the PAM stage will
respect the dendrogram in the sense an object can be PAM-assigned only to clusters that lie below it on
the branch that the object is merged into.
     See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{iteratePruningAndMerging}{Logical: should pruning of low-KME genes and module merging be iterated?
For backward compatibility, the default is \code{FALSE} but it setting it to \code{TRUE} may lead to
better-defined modules.}

  \item{minCoreKME}{ a number between 0 and 1. If a detected module does not have at least
\code{minModuleKMESize} genes with eigengene connectivity at least \code{minCoreKME}, the module is
disbanded (its genes are unlabeled and returned to the pool of genes waiting for mofule detection). }

  \item{minCoreKMESize}{ see \code{minCoreKME} above. }

  \item{minKMEtoStay}{ genes whose eigengene connectivity to their module eigengene is lower than
\code{minKMEtoStay} are removed from the module.}

  \item{impute}{ logical: should imputation be used for module eigengene calculation? See
\code{\link{moduleEigengenes}} for more details. }

 \item{trapErrors}{ logical: should errors in calculations be trapped? }
  \item{calibrateMergingSimilarities}{
Logical: should module eigengene similarities be calibrataed before calculating the consensus? Although
calibration is in principle desirable, the calibration methods currently available assume large data and do
not work very well on eigengene similarities.
}
  \item{mergeCutHeight}{
Dendrogram cut height for module merging. 
}
  \item{collectGarbage}{
Logical: should garbage be collected after some of the memory-intensive steps?
}
  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
  \item{\dots}{
Other arguments. Currently ignored.
}
}
\details{
This function calculates a consensus network with a flexible, possibly hierarchical consensus specification,
identifies (consensus) modules in the network, and calculates their eigengenes. "Blockwise" calculation is
available for large data sets for which a full network (TOM or adjacency matrix) would not fit into avilable
RAM. 

The input can be either several numerical data sets (expression etc) in the argument \code{multiExpr}
together with all necessary network construction options, or a pre-calculated network, typically the result
of a call to \code{\link{hierarchicalConsensusTOM}}. 

Steps in the network construction include the following: (1) optional 
filtering of variables (genes) and observations
(samples) that contain too many missing values or have zero variance; (2) optional pre-clustering to split
data into blocks of manageable size;
(3) calculation of adjacencies and optionally of TOMs in each individual data set;
(4) calculation of consensus network from the individual networks;
(5) hierarchical clustering and module identification;
(6) trimming of modules by removing genes with low correlation with the eigengene of the module; and 
(7) merging of modules whose eigengenes are strongly correlated.

Steps 1-4 (up to and including the calculation of consensus network from the individual networks) are
handled by the function \code{\link{hierarchicalConsensusTOM}}. 

Variables (genes) are clustered using average-linkage hierarchical clustering and modules are identified in the
resulting dendrogram by the Dynamic Hybrid tree cut. 

Found modules are trimmed of genes whose
consensus module membership kME (that is, correlation with module eigengene)
is less than \code{minKMEtoStay}.
Modules in which
fewer than \code{minCoreKMESize} genes have consensus KME higher than \code{minCoreKME}
are disbanded, i.e., their constituent genes are pronounced
unassigned.

After all blocks have been processed, the function checks whether there are genes whose KME in the module
they assigned is lower than KME to another module. If p-values of the higher correlations are smaller
than those of the native module by the factor \code{reassignThresholdPS} (in every set),
the gene is re-assigned to the closer module.

In the last step, modules whose eigengenes are highly correlated are merged. This is achieved by
clustering module eigengenes using the dissimilarity given by one minus their correlation,
cutting the dendrogram at the height \code{mergeCutHeight} and merging all modules on each branch. The
process is iterated until no modules are merged. See \code{\link{mergeCloseModules}} for more details on
module merging.

The module trimming and merging process is optionally iterated. Iterations are recommended but are (for now)
not the default for backward compatibility.

}
\value{
List with the following components:
\item{labels}{A numeric vector with one component per variable (gene), 
giving the module label of each variable (gene).  Label 0 is reserved for unassigned variables; module labels
are sequential and smaller numbers are used for larger modules.}

\item{unmergedLabels}{A numeric vector with one component per variable (gene),
giving the unmerged module label of each variable (gene), i.e., module labels before the call to module
merging.}

\item{colors}{A character vector with one component per variable (gene),
giving the module colors. The labels are mapped to colors using \code{\link{labels2colors}}.}

\item{unmergedColors}{A character vector with one component per variable (gene),
giving the unmerged module colors.}

\item{multiMEs}{Module eigengenes corresponding to the modules returned in \code{colors}, in multi-set
format. A vector of lists, one per set, containing eigengenes, proportion of variance explained and other
information. See \code{\link{multiSetMEs}} for a detailed description.}

\item{dendrograms}{A list with one component for each block of genes. Each component is the
hierarchical clustering dendrogram obtained by clustering the consensus gene dissimilarity in the
corresponding block. }

\item{consensusTOMInfo}{A list detailing various aspects of the consensus TOM. See
\code{\link{hierarchicalConsensusTOM}} for details.}

\item{blockInfo}{A list with information about blocks as well as the vriables and observations 
(genes and samples) retained after filtering out those with zero variance and too many missing values.}

\item{moduleIdentificationArguments}{A list with the module identification arguments supplied to this
function. Contains
\code{deepSplit},
\code{detectCutHeight},
\code{minModuleSize},
\code{maxCoreScatter},
\code{minGap},
\code{maxAbsCoreScatter},
\code{minAbsGap},
\code{minSplitHeight},
\code{useBranchEigennodeDissim},
\code{minBranchEigennodeDissim},
\code{minStabilityDissim},
\code{pamStage},
\code{pamRespectsDendro},
\code{minCoreKME},
\code{minCoreKMESize},
\code{minKMEtoStay},
\code{calibrateMergingSimilarities}, and
\code{mergeCutHeight}.}
}

\note{
If the input datasets have large numbers of genes, consider carefully the \code{maxBlockSize} as it
significantly affects the memory footprint (and whether the function will fail with a memory allocation
error). From a theoretical point of view it is advantageous to use blocks as large as possible; on the
other hand, using smaller blocks is substantially faster and often the only way to work with large
numbers of genes. As a rough guide, when 4GB of memory are available, blocks should be no larger than 8,000
genes; with 8GB one can handle some 13,000 genes; with 16GB around 20,000; and with 32GB around 30,000.
Depending on the operating system and its setup, these numbers may vary substantially.
}

\references{ 

Non-hierarchical consensus networks are described in Langfelder P, Horvath S (2007), Eigengene networks for
studying the relationships between co-expression modules. BMC Systems Biology 2007, 1:54. 

More in-depth discussion of selected topics can be found at http://www.peterlangfelder.com/ , and an FAQ at
https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html .
}


\author{
Peter Langfelder
}
\seealso{

\code{\link{hierarchicalConsensusTOM}} for calculation of hierarchical consensus networks (adjacency and
TOM), and a more detailed description of the calculation;

\code{\link[fastcluster]{hclust}} and \code{\link[dynamicTreeCut]{cutreeHybrid}} for hierarchical clustering
and the Dynamic Tree Cut branch cutting method;

\code{\link{mergeCloseModules}} for module merging;


\code{\link{blockwiseModules}} for an analogous analysis on a single data set.
}

\keyword{misc}% __ONLY ONE__ keyword per line
