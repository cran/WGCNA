\name{TOMsimilarityFromExpr}
\alias{TOMsimilarityFromExpr}
\title{ Topological overlap matrix }
\description{
  Calculation of the topological overlap matrix from given expression data. 
}
\usage{
TOMsimilarityFromExpr(
  datExpr, 
  weights = NULL,
  corType = "pearson", 
  networkType = "unsigned", 
  power = 6, 
  TOMType = "signed", 
  TOMDenom = "min",
  maxPOutliers = 1,
  quickCor = 0,
  pearsonFallback = "individual",
  cosineCorrelation = FALSE, 
  replaceMissingAdjacencies = FALSE,
  suppressTOMForZeroAdjacencies = FALSE,
  useInternalMatrixAlgebra = FALSE,
  nThreads = 0,
  verbose = 1, indent = 0)
}
\arguments{
  \item{datExpr}{ expression data. A data frame in which columns are genes and rows ar samples. NAs are
allowed, but not too many. }

  \item{weights}{optional observation weights for \code{datExpr} to be used in correlation calculation. 
A matrix of the same dimensions as
\code{datExpr}, containing non-negative weights.}

  \item{corType}{ character string specifying the correlation to be used. Allowed values are (unique
abbreviations of) \code{"pearson"} and \code{"bicor"}, corresponding to Pearson and bidweight
midcorrelation, respectively. Missing values are handled using the \code{pairwise.complete.obs} option. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{power}{ soft-thresholding power for netwoek construction. }


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

  \item{pearsonFallback}{Specifies whether the bicor calculation, if used, should revert to Pearson when median 
absolute deviation (mad) is zero. Recongnized values are (abbreviations of) 
\code{"none", "individual", "all"}. If set to
\code{"none"}, zero mad will result in \code{NA} for the corresponding correlation. 
If set to \code{"individual"}, Pearson calculation will be used only for columns that have zero mad. 
If set to \code{"all"}, the presence of a single zero mad will cause the whole variable to be treated in 
Pearson correlation manner (as if the corresponding \code{robust} option was set to \code{FALSE}). Has no
effect for Pearson correlation. See \code{\link{bicor}}.}

  \item{cosineCorrelation}{logical: should the cosine version of the correlation calculation be used? The 
cosine calculation differs from the standard one in that it does not subtract the mean. }

  \item{replaceMissingAdjacencies}{logical: should missing values in the calculation of adjacency be
replaced by 0?}

  \item{suppressTOMForZeroAdjacencies}{Logical: should TOM be set to zero for zero adjacencies?}

  \item{useInternalMatrixAlgebra}{Logical: should WGCNA's own, slow, matrix multiplication be used instead of
R-wide BLAS? Only useful for debugging.}

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
\value{
  A matrix holding the topological overlap.
}
\references{  
Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
}
\author{ Peter Langfelder}
\seealso{ \code{\link{TOMsimilarity}} }
\keyword{ misc }% __ONLY ONE__ keyword per line
