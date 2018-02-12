\name{sampledHierarchicalConsensusModules}
\alias{sampledHierarchicalConsensusModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical consensus module identification in sampled data
}
\description{
This function repeatedly resamples the samples (rows) in supplied data and identifies hierarchical consensus 
modules on the resampled data.
}
\usage{
sampledHierarchicalConsensusModules(
  multiExpr,
  multiWeights = NULL,

  networkOptions,
  consensusTree,

  nRuns,
  startRunIndex = 1,
  endRunIndex = startRunIndex + nRuns -1,
  replace = FALSE,
  fraction = if (replace) 1.0 else 0.63,
  randomSeed = 12345,
  checkSoftPower = TRUE,
  nPowerCheckSamples = 2000,
  individualTOMFilePattern = "individualTOM-Run.\%r-Set\%s-Block.\%b.RData",
  keepConsensusTOMs = FALSE,
  consensusTOMFilePattern = "consensusTOM-Run.\%r-\%a-Block.\%b.RData",
  skipUnsampledCalculation = FALSE,
  ...,
  verbose = 2, indent = 0,
  saveRunningResults = TRUE,
  runningResultsFile = "results.tmp.RData")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{
Expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes.
}

  \item{multiWeights}{ optional observation weights in the same format (and dimensions) as \code{multiExpr}.
These weights are used for correlation calculations with data in  \code{multiExpr}.}

  \item{networkOptions}{
A single list of class \code{\link{NetworkOptions}} giving options for network calculation for all of the
networks, or a \code{\link{multiData}} structure containing one such list for each input data set.
}
  \item{consensusTree}{
A list specifying the consensus calculation. See details.
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
  \item{individualTOMFilePattern}{Pattern for file names for files holding individual TOMs. The tags \code{"\%r, \%a,
\%b"} are replaced by run number, analysis name and block number, respectively. The TOM files are usually
temporary but can be retained, see \code{keepConsensusTOM} below.}

  \item{keepConsensusTOMs}{
Logical: should the (final) consensus TOMs of each sampled calculation be retained after the run ends? Note
that for large data sets (tens of thousands of nodes) the TOM files are rather large.
}
  \item{consensusTOMFilePattern}{

}
  \item{skipUnsampledCalculation}{
Logical: should a calculation on original (not resampled) data be skipped?
}
  \item{\dots}{
Other arguments to \code{\link{hierarchicalConsensusModules}}.
}

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

  \item{saveRunningResults}{
Logical: should the cumulative results be saved after each run on resampled data? 
}
  \item{runningResultsFile}{
File name of file in which to save running results into. In case of a parallel execution (say on several
nodes of a cluster), one should choose a unique name for each process to avoid overwriting the same file. 
}
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
\item{mods}{The output of the function \code{\link{hierarchicalConsensusModules}} on the resampled data.}
\item{samples}{Indices of the samples selected for the resampled data step for this run.}
\item{powers}{Actual soft-thresholding powers used in this run.}
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{hierarchicalConsensusModules}} for consensus networ analysis and module identification;

\code{\link{sampledBlockwiseModules}} for a similar resampling analysis for a single data set.
}
\keyword{misc}% __ONLY ONE__ keyword per line
