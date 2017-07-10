\name{blockwiseConsensusModules}
\alias{blockwiseConsensusModules}
\title{Find consensus modules across several datasets.}
\description{
  Perform network construction and consensus module detection across several datasets.
}
\usage{
blockwiseConsensusModules(
     multiExpr, 

     # Data checking options

     checkMissingData = TRUE,

     # Blocking options

     blocks = NULL, 
     maxBlockSize = 5000, 
     blockSizePenaltyPower = 5,
     nPreclusteringCenters = NULL,
     randomSeed = 12345,

     # TOM precalculation arguments, if available

     individualTOMInfo = NULL,
     useIndivTOMSubset = NULL,

     # Network construction arguments: correlation options

     corType = "pearson",
     maxPOutliers = 1,
     quickCor = 0,
     pearsonFallback = "individual", 
     cosineCorrelation = FALSE,

     # Adjacency function options

     power = 6, 
     networkType = "unsigned", 
     checkPower = TRUE,
     replaceMissingAdjacencies = FALSE,

     # Topological overlap options

     TOMType = "unsigned",
     TOMDenom = "min",

     # Save individual TOMs?

     saveIndividualTOMs = TRUE,
     individualTOMFileNames = "individualTOM-Set\%s-Block\%b.RData",

     # Consensus calculation options: network calibration

     networkCalibration = c("single quantile", "full quantile", "none"),

     # Simple quantile calibration options

     calibrationQuantile = 0.95,
     sampleForCalibration = TRUE, sampleForCalibrationFactor = 1000,
     getNetworkCalibrationSamples = FALSE,

     # Consensus definition

     consensusQuantile = 0,
     useMean = FALSE,
     setWeights = NULL,

     # Saving the consensus TOM

     saveConsensusTOMs = FALSE,
     consensusTOMFilePattern = "consensusTOM-block.\%b.RData",

     # Internal handling of TOMs

     useDiskCache = TRUE, chunkSize = NULL,
     cacheBase = ".blockConsModsCache",
     cacheDir = ".",

     # Alternative consensus TOM input from a previous calculation 

     consensusTOMInfo = NULL,

     # Basic tree cut options 

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
     minStabilityDissim = NULL,

     pamStage = TRUE,  pamRespectsDendro = TRUE,

     # Gene reassignment and trimming from a module, and module "significance" criteria

     reassignThresholdPS = 1e-4,
     trimmingConsensusQuantile = consensusQuantile,
     minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
     minKMEtoStay = 0.2,

     # Module eigengene calculation options

     impute = TRUE,
     trapErrors = FALSE,

     #Module merging options

     equalizeQuantilesForModuleMerging = FALSE,
     quantileSummaryForModuleMerging = "mean",
     mergeCutHeight = 0.15, 
     mergeConsensusQuantile = consensusQuantile,

     # Output options

     numericLabels = FALSE,

     # General options

     nThreads = 0,
     verbose = 2, indent = 0, ...)
}
\arguments{
  \item{multiExpr}{ expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

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
\code{nGenes/20} and \code{100*nGenes/maxBlockSize}, where \code{nGenes} is the nunber of genes (variables) in
\code{multiExpr}.}

 \item{randomSeed}{ integer to be used as seed for the random number generator before the function
starts. If a current seed exists, it is saved and restored upon exit. If \code{NULL} is given, the
function will not save and restore the seed. }

%%%%%%%%%%%%%%

  \item{individualTOMInfo}{ Optional data for TOM matrices in individual data sets. This object is returned by
the function \code{\link{blockwiseIndividualTOMs}}. If not given, appropriate topological overlaps will be
calculated using the network contruction options below. }

  \item{useIndivTOMSubset}{ If \code{individualTOMInfo} is given, this argument allows to only select a subset
of the individual set networks contained in \code{individualTOMInfo}. It should be a numeric vector giving the
indices of the individual sets to be used. Note that this argument is NOT applied to \code{multiExpr}. }


%%%%%%%%%%%%%%

  \item{corType}{ character string specifying the correlation to be used. Allowed values are (unique
abbreviations of) \code{"pearson"} and \code{"bicor"}, corresponding to Pearson and bidweight
midcorrelation, respectively. Missing values are handled using the \code{pariwise.complete.obs} option. }

  \item{maxPOutliers}{ only used for \code{corType=="bicor"}. Specifies the maximum percentile of data 
that can be considered outliers on either 
side of the median separately. For each side of the median, if
higher percentile than \code{maxPOutliers} is considered an outlier by the weight function based on
\code{9*mad(x)}, the width of the weight function is increased such that the percentile of outliers on
that side of the median equals \code{maxPOutliers}. Using \code{maxPOutliers=1} will effectively disable
all weight function broadening; using \code{maxPOutliers=0} will give results that are quite similar (but
not equal to) Pearson correlation. }


  \item{quickCor}{ real number between 0 and 1 that controls the handling of missing data in the
calculation of correlations. See details. }

  \item{pearsonFallback}{Specifies whether the bicor calculation, if used, should revert to Pearson when
median absolute deviation (mad) is zero. Recongnized values are (abbreviations of)
\code{"none", "individual", "all"}. If set to
\code{"none"}, zero mad will result in \code{NA} for the corresponding correlation.
If set to \code{"individual"}, Pearson calculation will be used only for columns that have zero mad.
If set to \code{"all"}, the presence of a single zero mad will cause the whole variable to be treated in
Pearson correlation manner (as if the corresponding \code{robust} option was set to \code{FALSE}). Has no
effect for Pearson correlation. See \code{\link{bicor}}.}

  \item{cosineCorrelation}{logical: should the cosine version of the correlation calculation be used? The 
cosine calculation differs from the standard one in that it does not subtract the mean. }

%%%%%%%%%%%%%%

  \item{power}{ soft-thresholding power for network construction. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{checkPower}{ logical: should basic sanity check be performed on the supplied \code{power}? If
you would like to experiment with unusual powers, set the argument to \code{FALSE} and proceed with
caution. }

  \item{replaceMissingAdjacencies}{logical: should missing values in the calculation of adjacency be replaced
by 0?}

  \item{TOMType}{ one of \code{"none"}, \code{"unsigned"}, \code{"signed"}. If \code{"none"}, adjacency
will be used for clustering. If \code{"unsigned"}, the standard TOM will be used (more generally, TOM
function will receive the adjacency as input). If \code{"signed"}, TOM will keep track of the sign of
correlations between neighbors. }

  \item{TOMDenom}{ a character string specifying the TOM variant to be used. Recognized values are 
\code{"min"} giving the standard TOM described in Zhang and Horvath (2005), and \code{"mean"} in which 
the \code{min} function in the denominator is replaced by \code{mean}. The \code{"mean"} may produce 
better results but at this time should be considered experimental.}
%The default mean denominator
%variant
%is preferrable and we recommend using it unless the user needs to reproduce older results obtained using
%the standard, minimum denominator TOM. }

%%%%%%%%%%%%%%%

  \item{saveIndividualTOMs}{logical: should individual TOMs be saved to disk for later use? }

  \item{individualTOMFileNames}{character string giving the file names to save individual TOMs into. The
following tags should be used to make the file names unique for each set and block: \code{\%s} will be
replaced by the set number; \code{\%N} will be replaced by the set name (taken from \code{names(multiExpr)})
if it exists, otherwise by set number; \code{\%b} will be replaced by the block number. If the file names turn
out to be non-unique, an error will be generated.}

%%%%%%%%%%%%%%

  \item{networkCalibration}{network calibration method. One of "single quantile", "full quantile", "none"
(or a unique abbreviation of one of them).}


%%%%%%%%%%%%%%

  \item{calibrationQuantile}{ if \code{networkCalibration} is \code{"single quantile"}, 
topological overlaps (or adjacencies if
TOMs are not computed) will be scaled such that their \code{calibrationQuantile} quantiles will agree. }

  \item{sampleForCalibration}{ if \code{TRUE}, calibration quantiles will be determined from a sample of network
similarities. Note that using all data can double the memory footprint of the function and the function
may fail. }

  \item{sampleForCalibrationFactor}{ determines the number of samples for calibration: the number is
 \code{1/calibrationQuantile * sampleForCalibrationFactor}. Should be set well above 1 to ensure accuracy of the
sampled quantile. }

 \item{getNetworkCalibrationSamples}{ logical: should samples used for TOM calibration be saved for future analysis?
This option is only available when \code{sampleForCalibration} is \code{TRUE}. }

%%%%%%%%%%%%%%

  \item{consensusQuantile}{ quantile at which consensus is to be defined. See details. }

  \item{useMean}{logical: should the consensus be determined from a (possibly weighted) mean across the
data sets rather than a quantile?}

  \item{setWeights}{Optional vector (one component per input set) of weights to be used for weighted mean
consensus. Only used when \code{useMean} above is \code{TRUE}.}

%%%%%%%%%%%%%%

  \item{saveConsensusTOMs}{ logical: should the consensus topological overlap matrices for each block be saved
and returned?  } 

  \item{consensusTOMFilePattern}{ character string containing the file namefiles containing the
consensus topological overlaps. The tag \code{\%b} will be replaced by the block number. If the resulting file
names are non-unique (for example, because the user gives a file name without a \code{\%b} tag), an error
will be generated.
These files are standard R data files and can be loaded using the \code{\link{load}}
function. } 

%%%%%%%%%%%%%%

  \item{useDiskCache}{ should calculated network similarities in individual sets be temporarilly saved
to disk? Saving to disk is somewhat slower than keeping all data in memory, but for large blocks and/or
many sets the memory footprint may be too big. }

  \item{chunkSize}{ network similarities are saved in smaller chunks of size \code{chunkSize}. }

  \item{cacheBase}{ character string containing the desired name for the cache files. The actual file
names will consists of \code{cacheBase} and a suffix to make the file names unique. }
  \item{cacheDir}{ character string containing the desired path for the cache files.}

%%%%%%%%%%%%%%

  \item{consensusTOMInfo}{optional list summarizing consensus TOM, output of \code{\link{consensusTOM}}. It
contains information about pre-calculated consensus TOM. Supplying this argument replaces TOM calculation,
so none of the individual or consensus TOM calculation arguments are taken into account.}

%%%%%%%%%%%%%%

  \item{deepSplit}{ integer value between 0 and 4. Provides a simplified control over how sensitive
module detection should be to module splitting, with 0 least and 4 most sensitive. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{detectCutHeight}{ dendrogram cut height for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  }

  \item{minModuleSize}{ minimum module size for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  }

  \item{checkMinModuleSize}{ logical: should sanity checks be performed on \code{minModuleSize}?}

%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%

  \item{reassignThresholdPS}{ per-set p-value ratio threshold for reassigning genes between modules. 
See Details. }

  \item{trimmingConsensusQuantile}{a number between 0 and 1 specifying the consensus quantile used for kME
calculation that determines module trimming according to the arguments below.}

  \item{minCoreKME}{ a number between 0 and 1. If a detected module does not have at least
\code{minModuleKMESize} genes with eigengene connectivity at least \code{minCoreKME}, the module is
disbanded (its genes are unlabeled and returned to the pool of genes waiting for mofule detection). }

  \item{minCoreKMESize}{ see \code{minCoreKME} above. }

  \item{minKMEtoStay}{ genes whose eigengene connectivity to their module eigengene is lower than
\code{minKMEtoStay} are removed from the module.}

%%%%%%%%%%%%%

  \item{impute}{ logical: should imputation be used for module eigengene calculation? See
\code{\link{moduleEigengenes}} for more details. }

 \item{trapErrors}{ logical: should errors in calculations be trapped? }

%%%%%%%%%%%%%

  \item{equalizeQuantilesForModuleMerging}{Logical: equalize quantiles of the module eigengene networks
before module merging? If \code{TRUE}, the quantiles of the eigengene correlation matrices (interpreted as a
single vectors of non-redundant components) will be equalized across the input data sets. Note that
although this seems like a reasonable option, it should be considered experimental and not necessarily
recommended.}

  \item{quantileSummaryForModuleMerging}{One of \code{"mean"} or \code{"median"}. 
If quantile equalization of the module eigengene networks is
performed, the resulting "normal" quantiles will be given by this function of the corresponding quantiles
across the input data sets.}


  \item{mergeCutHeight}{ dendrogram cut height for module merging. }

  \item{mergeConsensusQuantile}{consensus quantile for module merging. See \code{mergeCloseModules} for
details. }

  \item{numericLabels}{ logical: should the returned modules be labeled by colors (\code{FALSE}), or by
numbers (\code{TRUE})? }

%%%%%%%%%%%%%

  \item{nThreads}{ non-negative integer specifying the number of parallel threads to be used by certain
parts of correlation calculations. This option only has an effect on systems on which a POSIX thread
library is available (which currently includes Linux and Mac OSX, but excludes Windows).
If zero, the number of online processors will be used if it can be determined dynamically, otherwise
correlation calculations will use 2 threads. }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

  \item{...}{Other arguments. At present these can include \code{reproduceBranchEigennodeQuantileError} that
instructs the function to reproduce a bug in branch eigennode dissimilarity calculations for purposes if
reproducing old reults. }

}
\details{
  The function starts by optionally filtering out samples that have too many missing entries and genes
that have either too many missing entries or zero variance in at least one set. Genes that are filtered
out are left unassigned by the module detection. Returned eigengenes will contain \code{NA} in entries
corresponding to filtered-out samples. 

If \code{blocks} is not given and
the number of genes exceeds \code{maxBlockSize}, genes are pre-clustered into blocks using the function
\code{\link{consensusProjectiveKMeans}}; otherwise all genes are treated in a single block. 

For each block of genes, the network is constructed and (if requested) topological overlap is calculated
in each set. To minimize memory usage, calculated topological overlaps are optionally saved to disk in
chunks until they are needed again for the calculation of the consensus network topological overlap. 

Before calculation of the consensus Topological Overlap, individual TOMs are optionally calibrated.
Calibration methods include single quantile scaling and full quantile normalization. 

Single quantile
scaling raises individual TOM in sets 2,3,... to a power such that the quantiles given by 
\code{calibrationQuantile} agree with the quantile in set 1. Since the high TOMs are usually the most important
for module identification, the value of \code{calibrationQuantile} is close to (but not equal) 1. To speed up
quantile calculation, the quantiles can be determined on a randomly-chosen component subset of the TOM matrices.

Full quantile normalization, implemented in \code{\link[preprocessCore]{normalize.quantiles}}, adjusts the
TOM matrices such that all quantiles equal each other (and equal to the quantiles of the component-wise 
average of the individual TOM matrices). 

Note that network calibration is performed separately in each block, i.e., the normalizing transformation
may differ between blocks. This is necessary to avoid manipulating a full TOM in memory.

The consensus TOM is calculated as the component-wise \code{consensusQuantile} quantile of the individual
(set) TOMs; that is, for each gene pair (TOM entry), the \code{consensusQuantile} quantile across all input
sets. Alternatively, one can also use (weighted) component-wise mean across all imput data sets.
If requested, the consensus topological overlaps are saved to disk for later use.

Genes are then clustered using average linkage hierarchical clustering and modules are identified in the
resulting dendrogram by the Dynamic Hybrid tree cut. Found modules are trimmed of genes whose
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

 The argument \code{quick} specifies the precision of handling of missing data in the correlation 
calculations. Zero will cause all  
calculations to be executed precisely, which may be significantly slower than calculations without
missing data. Progressively higher values will speed up the
calculations but introduce progressively larger errors. Without missing data, all column means and
variances can be pre-calculated before the covariances are calculated. When missing data are present,
exact calculations require the column means and variances to be calculated for each covariance. The
approximate calculation uses the pre-calculated mean and variance and simply ignores missing data in the
covariance calculation. If the number of missing data is high, the pre-calculated means and variances may
be very different from the actual ones, thus potentially introducing large errors. 
The \code{quick} value times the
number of rows specifies the maximum difference in the
number of missing entries for mean and variance calculations on the one hand and covariance on the other
hand that will be tolerated before a recalculation is triggered. The hope is that if only a few missing
data are treated approximately, the error introduced will be small but the potential speedup can be
significant.
    
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

  \item{goodSamples}{ a list, with one component per input set. Each component is 
a logical vector with one entry per sample from the corresponding set. The entry indicates whether 
the sample in the set passed basic quality control criteria. }

  \item{goodGenes}{a logical vector with one entry per input gene indicating whether the gene passed
basic quality control criteria in all sets.}

  \item{dendrograms}{a list with one component for each block of genes. Each component is the
hierarchical clustering dendrogram obtained by clustering the consensus gene dissimilarity in the
corresponding block. }

  \item{TOMFiles}{ if \code{saveConsensusTOMs==TRUE}, 
a vector of character strings, one string per block, giving the file names of files
(relative to current directory) in which blockwise  topological overlaps were saved. }

  \item{blockGenes}{a list with one component for each block of genes. Each component is a vector giving
the indices (relative to the input \code{multiExpr}) of genes in the corresponding block. }

  \item{blocks}{if input \code{blocks} was given, its copy; otherwise a vector of length equal number of
genes giving the block label for each gene. Note that block labels are not necessarilly sorted in the
order in which the blocks were processed (since we do not require this for the input \code{blocks}). See
\code{blockOrder} below. }

  \item{blockOrder}{ a vector giving the order in which blocks were processed and in which
\code{blockGenes} above is returned. For example, \code{blockOrder[1]} contains the label of the
first-processed block. } 

  \item{originCount}{if the input \code{consensusQuantile==0}, this vector will contain counts of how
many times each set contributed the consensus gene similarity value. If the counts are highly unbalanced,
the consensus may be biased. }

  \item{networkCalibrationSamples}{if the input \code{getNetworkCalibrationSamples} is \code{TRUE}, this component is a
list with one component per block. Each component is again a list with two components:
\code{sampleIndex} contains indices of the distance structure in which TOM is stored that were sampled,
and \code{TOMSamples} is a matrix whose rows correspond to TOM samples and columns to individual set.
Hence, \code{networkCalibrationSamples[[blockNo]]$TOMSamples[index, setNo]} contains the TOM entry that
corresponds to element \code{networkCalibrationSamples[[blockNo]]$sampleIndex[index]} of the TOM distance
structure in block \code{blockNo} and set \code{setNo}. (For details on the distance structure, see
\code{\link{dist}}.)}

}

\note{ 
If the input datasets have large numbers of genes, consider carefully the \code{maxBlockSize} as it
significantly affects the memory footprint (and whether the function will fail with a memory allocation
error). From a theoretical point of view it is advantageous to use blocks as large as possible; on the
other hand, using smaller blocks is substantially faster and often the only way to work with large
numbers of genes. As a rough guide, it is unlikely a standard desktop
computer with 4GB memory or less will be able to work with blocks larger than 7000 genes. 

%Topological overlap calculations can be speeded up substantially (several 10-fold times on multi-core
%systems) if R is compiled with a dedicated BLAS (Basic Linear Algebra Subroutines) 
%library such as ATLAS or GotoBLAS and the package is compiled on your target system (which is always the
%case for Unix, Unix-like and Mac systems, but is normally not the case on Windows systems). 
}

\references{ Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between
co-expression modules. BMC Systems Biology 2007, 1:54 }

\author{Peter Langfelder}

\seealso{ 

\code{\link{goodSamplesGenesMS}} for basic quality control and filtering; 

\code{\link{adjacency}}, \code{\link{TOMsimilarity}} for network construction;

\code{\link{hclust}} for hierarchical clustering;

\code{\link[dynamicTreeCut]{cutreeDynamic}} for adaptive branch cutting in hierarchical clustering
dendrograms; 

\code{\link{mergeCloseModules}} for merging of close modules.
}

\keyword{ misc }% __ONLY ONE__ keyword per line
