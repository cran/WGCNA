\name{recutConsensusTrees}
\alias{recutConsensusTrees}

\title{ Repeat blockwise consensus module detection from pre-calculated data }
\description{
  Given consensus networks constructed for example using \code{\link{blockwiseConsensusModules}}, this
function (re-)detects modules in them by branch cutting of the corresponding dendrograms. If repeated
branch cuts of the same gene network dendrograms are desired, this function can save substantial time by
re-using already calculated networks and dendrograms.
}
\usage{
recutConsensusTrees(
  multiExpr,
  goodSamples, goodGenes,
  blocks,
  TOMFiles,
  dendrograms,
  corType = "pearson",
  networkType = "unsigned",
  deepSplit = 2,
  detectCutHeight = 0.995, minModuleSize = 20,
  checkMinModuleSize = TRUE,
  maxCoreScatter = NULL, minGap = NULL,
  maxAbsCoreScatter = NULL, minAbsGap = NULL,
  minSplitHeight = NULL, minAbsSplitHeight = NULL,

  useBranchEigennodeDissim = FALSE,
  minBranchEigennodeDissim = mergeCutHeight,

  pamStage = TRUE, pamRespectsDendro = TRUE,
  trimmingConsensusQuantile = 0,
  minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
  minKMEtoStay = 0.2,
  reassignThresholdPS = 1e-4,
  mergeCutHeight = 0.15,
  mergeConsensusQuantile = trimmingConsensusQuantile,
  impute = TRUE,
  trapErrors = FALSE,
  numericLabels = FALSE,
  verbose = 2, indent = 0)

}
\arguments{
  \item{multiExpr}{ expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

   \item{goodSamples}{ a list with one component per set. Each component is a logical vector specifying
which samples are considered "good" for the analysis. See \code{\link{goodSamplesGenesMS}}. }

   \item{goodGenes}{ a logical vector with length equal number of genes in \code{multiExpr} that
specifies which genes are considered "good" for the analysis. See \code{\link{goodSamplesGenesMS}}.  }

    \item{blocks}{ specification of blocks in which hierarchical clustering and module detection
should be performed. A numeric vector with one entry per gene
of \code{multiExpr} giving the number of the block to which the corresponding gene belongs. }

   \item{TOMFiles}{ a vector of character strings specifying file names in which the block-wise
topological overlaps are saved. }

   \item{dendrograms}{ a list of length equal the number of blocks, in which each component is a
hierarchical clustering dendrograms of the genes that belong to the block. }

  \item{corType}{ character string specifying the correlation to be used. Allowed values are (unique
abbreviations of) \code{"pearson"} and \code{"bicor"}, corresponding to Pearson and bidweight
midcorrelation, respectively. Missing values are handled using the \code{pariwise.complete.obs} option. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. Note that while no new networks are
computed in this function, this parameter affects the interpretation of correlations in this function. }

  \item{deepSplit}{ integer value between 0 and 4. Provides a simplified control over how sensitive
module detection should be to module splitting, with 0 least and 4 most sensitive. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{detectCutHeight}{ dendrogram cut height for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  }

  \item{minModuleSize}{ minimum module size for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  }

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

  \item{pamStage}{ logical.  If TRUE, the second (PAM-like) stage of module detection will be performed.
     See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{pamRespectsDendro}{Logical, only used when \code{pamStage} is \code{TRUE}.
If \code{TRUE}, the PAM stage will
respect the dendrogram in the sense an object can be PAM-assigned only to clusters that lie below it on
the branch that the object is merged into.
     See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{trimmingConsensusQuantile}{a number between 0 and 1 specifying the consensus quantile used for kME
calculation that determines module trimming according to the arguments below.}

  \item{minCoreKME}{ a number between 0 and 1. If a detected module does not have at least
\code{minModuleKMESize} genes with eigengene connectivity at least \code{minCoreKME}, the module is
disbanded (its genes are unlabeled and returned to the pool of genes waiting for mofule detection). }

  \item{minCoreKMESize}{ see \code{minCoreKME} above. }

  \item{minKMEtoStay}{ genes whose eigengene connectivity to their module eigengene is lower than
\code{minKMEtoStay} are removed from the module.}

  \item{reassignThresholdPS}{ per-set p-value ratio threshold for reassigning genes between modules. 
See Details. }

  \item{mergeCutHeight}{ dendrogram cut height for module merging. }

  \item{mergeConsensusQuantile}{consensus quantile for module merging. See \code{mergeCloseModules} for
details. }

  \item{impute}{ logical: should imputation be used for module eigengene calculation? See
\code{\link{moduleEigengenes}} for more details. }

 \item{trapErrors}{ logical: should errors in calculations be trapped? }

  \item{numericLabels}{ logical: should the returned modules be labeled by colors (\code{FALSE}), or by
numbers (\code{TRUE})? }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{

  For details on blockwise consensus module detection, see \code{\link{blockwiseConsensusModules}}. This
function implements the module detection subset of the functionality of 
\code{\link{blockwiseConsensusModules}}; network construction and clustering must be performed in
advance. The primary use of this function is to experiment with module detection settings without having
to re-execute long network and clustering calculations whose results are not affected by the cutting
parameters. 

This function takes as input the networks and dendrograms that are produced by
\code{\link{blockwiseConsensusModules}}.  Working block by block, 
modules are identified in the
dendrograms by the Dynamic Hybrid tree cut. 
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
    
}
\value{
  A list with the following components:

  \item{colors}{ module assignment of all input genes. A vector containing either character strings with
module colors (if input \code{numericLabels} was unset) or numeric module labels (if \code{numericLabels}
was set to \code{TRUE}). The color "grey" and the numeric label 0 are reserved for unassigned genes.  }

  \item{unmergedColors }{ module colors or numeric labels before the module merging step. }

  \item{multiMEs}{ module eigengenes corresponding to the modules returned in \code{colors}, in multi-set
format. A vector of lists, one per set, containing eigengenes, proportion of variance explained and other
information. See \code{\link{multiSetMEs}} for a detailed description. }

}
\note{ 
  Basic sanity checks are performed on given arguments, but it is left to the user's responsibility to
provide valid input.
}

\references{ Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between
co-expression modules. BMC Systems Biology 2007, 1:54 }

\author{Peter Langfelder}

\seealso{ 

\code{\link{blockwiseConsensusModules}} for the full blockwise modules calculation. Parts of its output
are natural input for this function.

\code{\link[dynamicTreeCut]{cutreeDynamic}} for adaptive branch cutting in hierarchical clustering
dendrograms; 

\code{\link{mergeCloseModules}} for merging of close modules.
}

\keyword{ misc }% __ONLY ONE__ keyword per line
