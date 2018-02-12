\name{blockwiseIndividualTOMs}
\alias{blockwiseIndividualTOMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of block-wise topological overlaps 
}
\description{
Calculates topological overlaps in the given (expression) data. If the number of variables (columns) in the
input data is too large, the data is first split using pre-clustering, then topological overlaps are
calculated in each block.
}
\usage{
blockwiseIndividualTOMs(
   multiExpr,
   multiWeights = NULL,

   # Data checking options

   checkMissingData = TRUE,

   # Blocking options

   blocks = NULL,
   maxBlockSize = 5000,
   blockSizePenaltyPower = 5,
   nPreclusteringCenters = NULL,
   randomSeed = 12345,

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
   suppressTOMForZeroAdjacencies = FALSE,

   # Topological overlap options

   TOMType = "unsigned",
   TOMDenom = "min",

   # Save individual TOMs? If not, they will be returned in the session.

   saveTOMs = TRUE,
   individualTOMFileNames = "individualTOM-Set\%s-Block\%b.RData",

   # General options

   nThreads = 0,
   useInternalMatrixAlgebra = FALSE,
   verbose = 2, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes.
}
  \item{multiWeights}{optional observation weights in the same format (and dimensions) as \code{multiExpr}.
These weights are used in correlation calculation.}

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

 \item{nPreclusteringCenters}{number of centers for pre-clustering. Larger numbers typically results in better
but slower pre-clustering. The default is \code{as.integer(min(nGenes/20, 100*nGenes/preferredSize))}
and is an attempt to arrive at a reasonable number given the resources available. }

 \item{randomSeed}{ integer to be used as seed for the random number generator before the function
starts. If a current seed exists, it is saved and restored upon exit. If \code{NULL} is given, the
function will not save and restore the seed. }

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

  \item{power}{ soft-thresholding power for netwoek construction. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{checkPower}{ logical: should basic sanity check be performed on the supplied \code{power}? If
you would like to experiment with unusual powers, set the argument to \code{FALSE} and proceed with
caution. }

  \item{replaceMissingAdjacencies}{logical: should missing values in calculated adjacency be replaced by 0?}

  \item{suppressTOMForZeroAdjacencies}{Logical: should TOM be set to zero for zero adjacencies?}

  \item{TOMType}{ one of \code{"none"}, \code{"unsigned"}, \code{"signed"}. If \code{"none"}, adjacency
will be used for clustering. If \code{"unsigned"}, the standard TOM will be used (more generally, TOM
function will receive the adjacency as input). If \code{"signed"}, TOM will keep track of the sign of
correlations between neighbors. Note that the \code{"unsigned"} vs. \code{"signed"} distinction is only
relevant when \code{networkType} is \code{"unsigned"}. When  \code{networkType} is \code{"signed"} or
\code{"signed hybrid"}, there is no difference between \code{TOMType="signed"} and \code{TOMType="unsigned".
} }

  \item{TOMDenom}{ a character string specifying the TOM variant to be used. Recognized values are
\code{"min"} giving the standard TOM described in Zhang and Horvath (2005), and \code{"mean"} in which
the \code{min} function in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
better results in certain special situations but at this time should be considered experimental.}
%The default mean denominator
%variant
%is preferrable and we recommend using it unless the user needs to reproduce older results obtained using
%the standard, minimum denominator TOM. }


  \item{saveTOMs}{logical: should calculated TOMs be saved to disk (\code{TRUE}) or returned in the return
value (\code{FALSE})? Returning calculated TOMs via the return value ay be more convenient bt not always
feasible if the matrices are too big to fit all in memory at the same time.
}

 \item{individualTOMFileNames}{character string giving the file names to save individual TOMs into. The
following tags should be used to make the file names unique for each set and block: \code{\%s} will be
replaced by the set number; \code{\%N} will be replaced by the set name (taken from \code{names(multiExpr)})
if it exists, otherwise by set number; \code{\%b} will be replaced by the block number. If the file names
turn out to be non-unique, an error will be generated.}

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

}
\details{
 The function starts by optionally filtering out samples that have too many missing entries and genes
that have either too many missing entries or zero variance in at least one set. Genes that are filtered
out are excluded from the TOM calculations. 

If \code{blocks} is not given and
the number of genes exceeds \code{maxBlockSize}, genes are pre-clustered into blocks using the function
\code{\link{consensusProjectiveKMeans}}; otherwise all genes are treated in a single block.

For each block of genes, the network is constructed and (if requested) topological overlap is calculated
in each set. The topological overlaps can be saved to disk as RData files, or returned directly within the
return value (see below). Note that the matrices can be big and returning them within the return value can
quickly exhaust the system's memory. In particular, if the block-wise calculation is necessary, it is
nearly certain that returning all matrices via the return value will be impossible.
}

\value{
A list with the following components:

  \item{actualTOMFileNames}{Only returned if input \code{saveTOMs} is \code{TRUE}. A matrix of character
strings giving the file names in which each block TOM is saved. Rows correspond to data sets and columns to
blocks.}

  \item{TOMSimilarities}{Only returned if input \code{saveTOMs} is \code{FALSE}. A list in which each
component corresponds to one block. Each component is a matrix of dimensions (N times (number of sets)), where
N is the length of a distance structure corresponding to the block. That is, if the block contains n genes,
N=n*(n-1)/2. Each column of the matrix contains the topological overlap of variables in the corresponding set (
and the corresponding block), arranged as a distance structure. Do note however that the topological overlap
is a similarity (not a distance). }

  \item{blocks}{if input \code{blocks} was given, its copy; otherwise a vector of length equal number of
genes giving the block label for each gene. Note that block labels are not necessarilly sorted in the
order in which the blocks were processed (since we do not require this for the input \code{blocks}). See
\code{blockOrder} below. }

  \item{blockGenes}{a list with one component for each block of genes. Each component is a vector giving
the indices (relative to the input \code{multiExpr}) of genes in the corresponding block. }

  \item{goodSamplesAndGenes}{if input 
\code{checkMissingData} is \code{TRUE}, the output of the function \code{\link{goodSamplesGenesMS}}. 
A list with components
\code{goodGenes} (logical vector indicating which genes passed the missing data filters), \code{goodSamples}
(a list of logical vectors indicating which samples passed the missing data filters in each set), and
\code{allOK} (a logical indicating whether all genes and all samples passed the filters). See
\code{\link{goodSamplesGenesMS}} for more details. If \code{checkMissingData} is \code{FALSE},
\code{goodSamplesAndGenes} contains a list of the same type but indicating that all genes and all samples
passed the missing data filters.} 

The following components are present mostly to streamline the interaction of this function with
\code{\link{blockwiseConsensusModules}}.

  \item{nGGenes}{ Number of genes that passed missing data filters (if input  
\code{checkMissingData} is \code{TRUE}), or the number of all genes (if \code{checkMissingData} is
\code{FALSE}).}

  \item{gBlocks}{ the vector \code{blocks} (above), restricted to good genes only. }

  \item{nThreads}{ number of threads used to calculate correlation and TOM matrices. }

  \item{saveTOMs}{ logical: were calculated matrices saved in files (\code{TRUE}) or returned in the
return value (\code{FALSE})?}

  \item{intNetworkType, intCorType}{integer codes for network and  correlation type. }

  \item{nSets}{number of sets in input data.}

  \item{setNames}{the \code{names} attribute of input \code{multiExpr}.}
}

\references{
For a general discussion of the weighted network formalism, see

Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network Analysis",
Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17

The blockwise approach is briefly described in the article describing this package,

Langfelder P, Horvath S (2008) "WGCNA: an R package for weighted correlation network analysis". 
BMC Bioinformatics 2008, 9:559

}
\author{
Peter Langfelder
}

\seealso{

\code{\link{blockwiseConsensusModules}}

}
\keyword{misc}
