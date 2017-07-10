\name{sampledBlockwiseModules}
\alias{sampledBlockwiseModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Blockwise module identification in sampled data
}
\description{
This function repeatedly resamples the samples (rows) in supplied data and identifies 
modules on the resampled data.
}
\usage{
sampledBlockwiseModules(
  datExpr,
  nRuns,
  startRunIndex = 1,
  endRunIndex = startRunIndex + nRuns -1,
  replace = FALSE,
  fraction = if (replace) 1.0 else 0.63,
  randomSeed = 12345,
  checkSoftPower = TRUE,
  nPowerCheckSamples = 2000,
  skipUnsampledCalculation = FALSE,
  corType = "pearson",
  power = 6,
  networkType = "unsigned",
  saveTOMs = FALSE,
  saveTOMFileBase = "TOM",
  ...,
  verbose = 2, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{
 Expression data. A matrix (preferred) or 
data frame in which columns are genes and rows ar samples. 
}
  \item{nRuns}{
Number of network construction and module identification runs.
}
  \item{startRunIndex}{
Number to be assigned to the start run. The run number or index is used to make saved files unique; it has no
effect on the actual results of the run.
}
  \item{endRunIndex}{
Number (index) of the last run. If given, \code{nRuns} is ignored.
}
  \item{replace}{
Logical: should samples (observations or rows in entries in \code{multiExpr}) be sampled with replacement?
}
  \item{fraction}{
Fraction of samples to sample for each run.
}
  \item{randomSeed}{
Integer specifying the random seed. If non-NULL, the random number generator state is saved before the seed is set
and restored at the end of the function. If \code{NULL}, the random number generator state is not changed nor
saved at the start, and not restored at the end.
}
  \item{checkSoftPower}{
Logical: should the soft-tresholding power be adjusted to approximately match the connectivity distribution
of the sampled data set and the full data set?
}
  \item{nPowerCheckSamples}{
Number of genes to be sampled from the full data set to calculate connectivity and match soft-tresholding
powers.
}
  \item{skipUnsampledCalculation}{
Logical: should a calculation on original (not resampled) data be skipped?
}
 \item{corType}{Character string specifying the correlation to be used. Allowed values are (unique
abbreviations of) \code{"pearson"} and \code{"bicor"}, corresponding to Pearson and bidweight
midcorrelation, respectively. Missing values are handled using the \code{pairwise.complete.obs} option. }
  \item{power}{ Soft-thresholding power for network construction. }

 \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{saveTOMs}{
Logical: should the networks (topological overlaps) be saved for each run? Note
that for large data sets (tens of thousands of nodes) the TOM files are rather large.
}
  \item{saveTOMFileBase}{
Character string giving the base of the file names for TOMs. The actual file names will consist of a
concatenation of \code{saveTOMFileBase} and \code{"-run-<run
number>-Block-<block number>.RData"}. 
}
  \item{\dots}{
Other arguments to \code{\link{blockwiseModules}}.
}

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
For each run, samples (but not genes) are randomly sampled to obtain a perturbed data set; a full network
analysis and module identification is carried out, and the results are returned in a list with one component
per run.

For each run, the soft-thresholding power can optionally be adjusted such that the mean adjacency in the
re-sampled data set equals the mean adjacency in the original data. 
}
\value{
A list with one component per run. Each component is a list with the following components:
\item{mods}{The output of the function \code{\link{blockwiseModules}} applied to a resampled data set.}
\item{samples}{Indices of the samples selected for the resampled data step for this run.}
\item{powers}{Actual soft-thresholding powers used in this run.}
}
\references{
An application of this function is described in the motivational example section of 

Langfelder P, Horvath S (2012) Fast R Functions for Robust Correlations and Hierarchical Clustering.
Journal of Statistical Software 46(11) 1-17; PMID: 23050260 PMCID: PMC3465711 
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{blockwiseModules}} for the underlying network analysis and module identification;

\code{\link{sampledHierarchicalConsensusModules}} for a similar resampling analysis of consensus networks.

}
\keyword{misc}% __ONLY ONE__ keyword per line
