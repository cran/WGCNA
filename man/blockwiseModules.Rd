\name{blockwiseModules}
\alias{blockwiseModules}
\title{ Automatic network construction and module detection }
\description{
 This function performs automatic network construction and module detection on large expression datasets
in a block-wise manner.  
}
\usage{
blockwiseModules(
  # Input data

  datExpr, 
  weights = NULL,

  # Data checking options

  checkMissingData = TRUE,

  # Options for splitting data into blocks

  blocks = NULL,
  maxBlockSize = 5000,
  blockSizePenaltyPower = 5,
  nPreclusteringCenters = as.integer(min(ncol(datExpr)/20, 
                             100*ncol(datExpr)/maxBlockSize)),
  randomSeed = 12345,

 # load TOM from previously saved file?

  loadTOM = FALSE,

  # Network construction arguments: correlation options

  corType = "pearson",
  maxPOutliers = 1, 
  quickCor = 0,
  pearsonFallback = "individual",
  cosineCorrelation = FALSE,

  # Adjacency function options

  power = 6,
  networkType = "unsigned",
  replaceMissingAdjacencies = FALSE,
  suppressTOMForZeroAdjacencies = FALSE,

  # Topological overlap options

  TOMType = "signed",
  TOMDenom = "min",

  # Saving or returning TOM

  getTOMs = NULL,
  saveTOMs = FALSE, 
  saveTOMFileBase = "blockwiseTOM",

  # Basic tree cut options

  deepSplit = 2,
  detectCutHeight = 0.995, 
  minModuleSize = min(20, ncol(datExpr)/2 ),

  # Advanced tree cut options

  maxCoreScatter = NULL, minGap = NULL,
  maxAbsCoreScatter = NULL, minAbsGap = NULL,
  minSplitHeight = NULL, minAbsSplitHeight = NULL,

  useBranchEigennodeDissim = FALSE,
  minBranchEigennodeDissim = mergeCutHeight,

  stabilityLabels = NULL,
  stabilityCriterion = c("Individual fraction", "Common fraction"),
  minStabilityDissim = NULL,

  pamStage = TRUE, pamRespectsDendro = TRUE,

  # Gene reassignment, module trimming, and module "significance" criteria

  reassignThreshold = 1e-6,
  minCoreKME = 0.5, 
  minCoreKMESize = minModuleSize/3,
  minKMEtoStay = 0.3,

  # Module merging options

  mergeCutHeight = 0.15, 
  impute = TRUE, 
  trapErrors = FALSE, 

  # Output options

  numericLabels = FALSE,

  # Options controlling behaviour

  nThreads = 0,
  useInternalMatrixAlgebra = FALSE,
  verbose = 0, indent = 0,
  ...)

}
\arguments{
  \item{datExpr}{ Expression data. A matrix (preferred) or 
data frame in which columns are genes and rows ar samples. NAs are
allowed, but not too many. See \code{checkMissingData} below and details.} 

  \item{weights}{optional observation weights in the same format (and dimensions) as \code{datExpr}.
These weights are used in correlation calculation.}

  \item{checkMissingData}{logical: should data be checked for excessive numbers of missing entries in
genes and samples, and for genes with zero variance? See details. }

  \item{blocks}{optional specification of blocks in which hierarchical clustering and module detection
should be performed. If given, must be a numeric vector with one entry per column (gene) 
of \code{exprData} giving the number of the block to which the corresponding gene belongs. }

  \item{maxBlockSize}{integer giving maximum block size for module detection. Ignored if \code{blocks}
above is non-NULL. Otherwise, if the number of genes in \code{datExpr} exceeds \code{maxBlockSize}, genes
will be pre-clustered into blocks whose size should not exceed \code{maxBlockSize}. }

 \item{ blockSizePenaltyPower}{number specifying how strongly blocks should be penalized for exceeding the
maximum size. Set to a lrge number or \code{Inf} if not exceeding maximum block size is very important.}

 \item{nPreclusteringCenters}{number of centers for pre-clustering. Larger numbers typically results in better
but slower pre-clustering. }
 
 \item{randomSeed}{integer to be used as seed for the random number generator before the function
starts. If a current seed exists, it is saved and restored upon exit. If \code{NULL} is given, the 
function will not save and restore the seed. }

  \item{loadTOM}{logical: should Topological Overlap Matrices be loaded from previously saved files (\code{TRUE}) 
or calculated (\code{FALSE})? It may be useful to load previously saved TOM matrices if these have been
calculated previously, since TOM calculation is often the most computationally expensive part of network
construction and module identification. See \code{saveTOMs} and \code{saveTOMFileBase} below for when and how TOM
files are saved, and what the file names are. If \code{loadTOM} is \code{TRUE} but the files cannot be
found, or do not contain the correct TOM data, TOM will be recalculated.}

  \item{corType}{character string specifying the correlation to be used. Allowed values are (unique
abbreviations of) \code{"pearson"} and \code{"bicor"}, corresponding to Pearson and bidweight
midcorrelation, respectively. Missing values are handled using the \code{pairwise.complete.obs} option. }

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
effect for Pearson correlation.  See \code{\link{bicor}}.}

  \item{cosineCorrelation}{logical: should the cosine version of the correlation calculation be used? The
cosine calculation differs from the standard one in that it does not subtract the mean. } 
  
  \item{power}{ soft-thresholding power for network construction. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{replaceMissingAdjacencies}{logical: should missing values in the calculation of adjacency be replaced
by 0?}

  \item{suppressTOMForZeroAdjacencies}{Logical: should TOM be set to zero for zero adjacencies?}

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

  \item{getTOMs}{ deprecated, please use saveTOMs below. }

  \item{saveTOMs}{ logical: should the consensus topological overlap matrices for each block be saved
and returned?  } 

  \item{saveTOMFileBase}{ character string containing the file name base for files containing the
consensus topological overlaps. The full file names have \code{"block.1.RData"}, \code{"block.2.RData"}
etc. appended. These files are standard R data files and can be loaded using the \code{\link{load}}
function. }


  \item{deepSplit}{ integer value between 0 and 4. Provides a simplified control over how sensitive
module detection should be to module splitting, with 0 least and 4 most sensitive. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

  \item{detectCutHeight}{ dendrogram cut height for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  }

  \item{minModuleSize}{ minimum module size for module detection. See
\code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.  }

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

  \item{minCoreKME}{ a number between 0 and 1. If a detected module does not have at least
\code{minModuleKMESize} genes with eigengene connectivity at least \code{minCoreKME}, the module is
disbanded (its genes are unlabeled and returned to the pool of genes waiting for mofule detection). }

  \item{minCoreKMESize}{ see \code{minCoreKME} above. }

  \item{minKMEtoStay}{ genes whose eigengene connectivity to their module eigengene is lower than
\code{minKMEtoStay} are removed from the module.}

  \item{reassignThreshold}{ p-value ratio threshold for reassigning genes between modules. See Details. }

  \item{mergeCutHeight}{ dendrogram cut height for module merging. }

  \item{impute}{ logical: should imputation be used for module eigengene calculation? See
\code{\link{moduleEigengenes}} for more details. }

  \item{trapErrors}{ logical: should errors in calculations be trapped? }

  \item{numericLabels}{ logical: should the returned modules be labeled by colors (\code{FALSE}), or by
numbers (\code{TRUE})? }

  \item{nThreads}{ non-negative integer specifying the number of parallel threads to be used by certain
parts of correlation calculations. This option only has an effect on systems on which a POSIX thread
library is available (which currently includes Linux and Mac OSX, but excludes Windows). 
If zero, the number of online processors will be used if it can be determined dynamically, otherwise
correlation calculations will use 2 threads. }

  \item{useInternalMatrixAlgebra}{Logical: should WGCNA's own, slow, matrix multiplication be used instead of
R-wide BLAS? Only useful for debugging.}  

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

  \item{...}{Other arguments.}
}
\details{
Before module detection starts, genes and samples are optionally checked for the presence of \code{NA}s.
Genes and/or samples that have too many \code{NA}s are flagged as bad and removed from the analysis; bad
genes will be automatically labeled as unassigned, while the returned eigengenes will have \code{NA}
entries for all bad samples. 

If \code{blocks} is not given and
the number of genes exceeds \code{maxBlockSize}, genes are pre-clustered into blocks using the function
\code{\link{projectiveKMeans}}; otherwise all genes are treated in a single block.

For each block of genes, the network is constructed and (if requested) topological overlap is calculated.
If requested, the topological overlaps are returned as part of the return value list.
Genes are then clustered using average linkage hierarchical clustering and modules are identified in the
resulting dendrogram by the Dynamic Hybrid tree cut. Found modules are trimmed of genes whose
correlation with module eigengene (KME) is less than \code{minKMEtoStay}. Modules in which
fewer than \code{minCoreKMESize} genes have KME higher than \code{minCoreKME} 
are disbanded, i.e., their constituent genes are pronounced
unassigned. 

After all blocks have been processed, the function checks whether there are genes whose KME in the module
they assigned is lower than KME to another module. If p-values of the higher correlations are smaller
than those of the native module by the factor \code{reassignThresholdPS},
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

  \item{colors }{ a vector of color or numeric module labels for all genes.}

  \item{unmergedColors }{ a vector of color or numeric module labels for all genes before module merging.}

  \item{MEs }{ a data frame containing module eigengenes of the found modules (given by \code{colors}).}

  \item{goodSamples}{numeric vector giving indices of good samples, that is samples that do not have too
many missing entries. }

  \item{goodGenes}{ numeric vector giving indices of good genes, that is genes that do not have too
many missing entries.}

  \item{dendrograms}{ a list whose components conatain hierarchical clustering dendrograms of genes 
in each block. }

  \item{TOMFiles}{ if \code{saveTOMs==TRUE},
a vector of character strings, one string per block, giving the file names of files
(relative to current directory) in which blockwise  topological overlaps were saved. }

  \item{blockGenes}{ a list whose components give the indices of genes in each block. }

  \item{blocks}{if input \code{blocks} was given, its copy; otherwise a vector of length equal number of
genes giving the block label for each gene. Note that block labels are not necessarilly sorted in the
order in which the blocks were processed (since we do not require this for the input \code{blocks}). See
\code{blockOrder} below. }

  \item{blockOrder}{ a vector giving the order in which blocks were processed and in which
\code{blockGenes} above is returned. For example, \code{blockOrder[1]} contains the label of the
first-processed block. }

  \item{MEsOK}{logical indicating whether the module eigengenes were calculated without errors. }
}

\note{
significantly affects the memory footprint (and whether the function will fail with a memory allocation
error). From a theoretical point of view it is advantageous to use blocks as large as possible; on the
other hand, using smaller blocks is substantially faster and often the only way to work with large
numbers of genes. As a rough guide, it is unlikely a standard desktop
computer with 4GB memory or less will be able to work with blocks larger than 8000 genes.
}

\references{Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 }

\author{ Peter Langfelder}
% R documentation directory.
\seealso{ 

\code{\link{goodSamplesGenes}} for basic quality control and filtering;

\code{\link{adjacency}}, \code{\link{TOMsimilarity}} for network construction;

\code{\link[stats]{hclust}} for hierarchical clustering;
  
\code{\link[dynamicTreeCut]{cutreeDynamic}} for adaptive branch cutting in hierarchical clustering
dendrograms;

\code{\link{mergeCloseModules}} for merging of close modules.

}

\keyword{ misc }
