\name{consensusTOM}
\alias{consensusTOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Consensus network (topological overlap).}
\description{
Calculation of a consensus network (topological overlap).
}
\usage{
consensusTOM(
      # Supply either ...
      # ... information needed to calculate individual TOMs
      multiExpr,

      # Data checking options
      checkMissingData = TRUE,

      # Blocking options
      blocks = NULL,
      maxBlockSize = 5000,
      blockSizePenaltyPower = 5,
      nPreclusteringCenters = NULL,
      randomSeed = 54321,

      # Network construction arguments: correlation options

      corType = "pearson",
      maxPOutliers = 1,
      quickCor = 0,
      pearsonFallback = "individual",
      cosineCorrelation = FALSE,
      replaceMissingAdjacencies = FALSE,

      # Adjacency function options

      power = 6,
      networkType = "unsigned",
      checkPower = TRUE,

      # Topological overlap options

      TOMType = "unsigned",
      TOMDenom = "min",
      suppressNegativeTOM = FALSE,

      # Save individual TOMs?

      saveIndividualTOMs = TRUE,
      individualTOMFileNames = "individualTOM-Set\%s-Block\%b.RData",

      # ... or individual TOM information

      individualTOMInfo = NULL,
      useIndivTOMSubset = NULL,

   ##### Consensus calculation options 

      useBlocks = NULL,

      networkCalibration = c("single quantile", "full quantile", "none"),

      # Save calibrated TOMs?
      saveCalibratedIndividualTOMs = FALSE,
      calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set\%s-Block\%b.RData",

      # Simple quantile calibration options
      calibrationQuantile = 0.95,
      sampleForCalibration = TRUE, sampleForCalibrationFactor = 1000,
      getNetworkCalibrationSamples = FALSE,

      # Consensus definition
      consensusQuantile = 0,
      useMean = FALSE,
      setWeights = NULL,

      # Return options
      saveConsensusTOMs = TRUE,
      consensusTOMFilePattern = "consensusTOM-Block\%b.RData",
      returnTOMs = FALSE,

      # Internal handling of TOMs
      useDiskCache = NULL, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      nThreads = 1,

      # Diagnostic messages
      verbose = 1,
      indent = 0)
}
%- maybe also 'usage' for other objects documented here.
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

 \item{ blockSizePenaltyPower}{number specifying how strongly blocks should be penalized for exceeding the
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

 \item{power}{ soft-thresholding power for network construction. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{checkPower}{ logical: should basic sanity check be performed on the supplied \code{power}? If
you would like to experiment with unusual powers, set the argument to \code{FALSE} and proceed with
caution. }

  \item{replaceMissingAdjacencies}{logical: should missing values in the calculation of adjacency be replaced
by 0?}

  \item{TOMType}{ one of \code{"none"}, \code{"unsigned"}, \code{"signed"}, \code{"signed Nowick"},
\code{"unsigned 2"}, \code{"signed 2"} and \code{"signed Nowick 2"}. If \code{"none"}, adjacency
will be used for clustering. See \code{\link{TOMsimilarityFromExpr}} for details.}

  \item{TOMDenom}{ a character string specifying the TOM variant to be used. Recognized values are
\code{"min"} giving the standard TOM described in Zhang and Horvath (2005), and \code{"mean"} in which
the \code{min} function in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
better results but at this time should be considered experimental.}

 \item{suppressNegativeTOM}{Logical: should the result be set to zero when negative? Negative TOM values can occur when
\code{TOMType} is \code{"signed Nowick"}.}


%%%%%%%%%%%%%%%

  \item{saveIndividualTOMs}{logical: should individual TOMs be saved to disk for later use? }

  \item{individualTOMFileNames}{character string giving the file names to save individual TOMs into. The
following tags should be used to make the file names unique for each set and block: \code{\%s} will be
replaced by the set number; \code{\%N} will be replaced by the set name (taken from \code{names(multiExpr)})
if it exists, otherwise by set number; \code{\%b} will be replaced by the block number. If the file names turn
out to be non-unique, an error will be generated.}

%%%%%%%%%%%%%%

  \item{individualTOMInfo}{ Optional data for TOM matrices in individual data sets. This object is returned by
the function \code{\link{blockwiseIndividualTOMs}}. If not given, appropriate topological overlaps will be
calculated using the network contruction options below. }

  \item{useIndivTOMSubset}{ If \code{individualTOMInfo} is given, this argument allows to only select a subset
of the individual set networks contained in \code{individualTOMInfo}. It should be a numeric vector giving the
indices of the individual sets to be used. Note that this argument is NOT applied to \code{multiExpr}. }

  \item{useBlocks}{optional specification of blocks that should be used for the calcualtions. The default is to
use all blocks.

}
%%%%%%%%%%%%%%

  \item{networkCalibration}{network calibration method. One of "single quantile", "full quantile", "none"
(or a unique abbreviation of one of them).}

  \item{saveCalibratedIndividualTOMs}{logical: should the calibrated individual TOMs be saved?
}
  \item{calibratedIndividualTOMFilePattern}{pattern of file names for saving calibrated individual TOMs.}

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

  \item{getNetworkCalibrationSamples}{logical: should the sampled values used for network calibration be
returned?}

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
  \item{returnTOMs}{logical: should calculated consensus TOM(s) be returned?
}

  \item{useDiskCache}{ should calculated network similarities in individual sets be temporarilly saved
to disk? Saving to disk is somewhat slower than keeping all data in memory, but for large blocks and/or
many sets the memory footprint may be too big. If not given (the default), the function will determine
the need of caching based on the size of the data. See \code{chunkSize} below for additional information. }

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

  \item{nThreads}{ non-negative integer specifying the number of parallel threads to be used by certain
parts of correlation calculations. This option only has an effect on systems on which a POSIX thread
library is available (which currently includes Linux and Mac OSX, but excludes Windows).
If zero, the number of online processors will be used if it can be determined dynamically, otherwise
correlation calculations will use 2 threads. }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

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
\code{calibrationQuantile} agree with the quantile in set 1. Since the high TOMs are usually the most
important
for module identification, the value of \code{calibrationQuantile} is close to (but not equal) 1. To speed up
quantile calculation, the quantiles can be determined on a randomly-chosen component subset of the TOM
matrices.

Full quantile normalization, implemented in \code{\link[preprocessCore]{normalize.quantiles}}, adjusts the
TOM matrices such that all quantiles equal each other (and equal to the quantiles of the component-wise
average of the individual TOM matrices).

Note that network calibration is performed separately in each block, i.e., the normalizing transformation
may differ between blocks. This is necessary to avoid manipulating a full TOM in memory.

The consensus TOM is calculated as the component-wise \code{consensusQuantile} quantile of the individual
(set) TOMs; that is, for each gene pair (TOM entry), the \code{consensusQuantile} quantile across all input
sets. Alternatively, one can also use (weighted) component-wise mean across all imput data sets.
If requested, the consensus topological overlaps are saved to disk for later use.
}
\value{
List with the following components:

\item{consensusTOM}{only present if input \code{returnTOMs} is \code{TRUE}. A list containing consensus TOM
for each block, stored as a distance structure.}

\item{TOMFiles}{only present if input \code{saveConsensusTOMs} is \code{TRUE}. A vector of file names, one for
each block, in which the TOM for the corresponding block is stored. TOM is saved as a distance structure to
save space.}

\item{saveConsensusTOMs}{a copy of the input{saveConsensusTOMs}.}

\item{individualTOMInfo}{information about individual set TOMs. A copy of the input \code{individualTOMInfo}
if given; otherwise the result of calling \code{blockwiseIndividualTOMs}. See \code{blockwiseIndividualTOMs} for
details.}

Further components are retained for debugging and/or convenience.

\item{useIndivTOMSubset}{a copy of the input \code{useIndivTOMSubset}.}

\item{goodSamplesAndGenes}{a list containing information about which samples and genes are "good" in the sense
that they do not contain more than a certain fraction of missing data and (for genes) have non-zero variance.
See \code{\link{goodSamplesGenesMS}} for details.}

\item{nGGenes}{number of "good" genes in \code{goodSamplesGenes} above. }

\item{nSets}{number of input sets.}

\item{saveCalibratedIndividualTOMs}{a copy of the input \code{saveCalibratedIndividualTOMs}.}

\item{calibratedIndividualTOMFileNames}{if input \code{saveCalibratedIndividualTOMs} is \code{TRUE}, this
component will contain the file names of calibrated individual networks. The file names are arranged in a
character matrix with each row corresponding to one input set and each column to one block.}

\item{networkCalibrationSamples}{if input \code{getNetworkCalibrationSamples} is \code{TRUE}, a list with one
component per block. Each component is in turn a list with two components: \code{sampleIndex} is a vector
contain the indices of the TOM samples (the indices refer to a flattened distance structure), and
\code{TOMSamples} is a matrix of TOM samples with each row corresponding to a sample in \code{sampleIndex},
and each column to one input set.}

\item{consensusQuantile}{a copy of the input \code{consensusQuantile}.}

\item{originCount}{A vector of length \code{nSets} that
contains, for each set, the number of (calibrated) elements that were less than or equal the consensus for that element.}

}
\references{
WGCNA methodology has been described in

Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network Analysis",
Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834

The original reference for the WGCNA package is 

Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC
Bioinformatics 2008, 9:559 PMID: 19114008

For consensus modules, see 

Langfelder P, Horvath S (2007) "Eigengene networks for studying the relationships between
co-expression modules", BMC Systems Biology 2007, 1:54

This function uses quantile normalization described, for example, in

Bolstad BM1, Irizarry RA, Astrand M, Speed TP (2003) "A comparison of normalization methods for high density
oligonucleotide array data based on variance and bias", Bioinformatics. 2003 Jan 22;19(2):1 
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{blockwiseIndividualTOMs}}  for calculation of topological overlaps across multiple sets.
}
\keyword{misc}% __ONLY ONE__ keyword per line
