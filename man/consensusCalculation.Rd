\name{consensusCalculation}
\alias{consensusCalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of a (single) consenus with optional data calibration.
}
\description{
This function calculates a single consensus from given individual data, optionally first calibrating the
individual data to make them comparable.
}
\usage{
consensusCalculation(
  individualData,
  consensusOptions,

  useBlocks = NULL,
  randomSeed = NULL,
  saveCalibratedIndividualData = FALSE,
  calibratedIndividualDataFilePattern = "calibratedIndividualData-\%a-Set\%s-Block\%b.RData",

  # Return options: the data can be either saved or returned but not both.
  saveConsensusData = NULL,
  consensusDataFileNames = "consensusData-\%a-Block\%b.RData",
  getCalibrationSamples= FALSE,

  # Internal handling of data
  useDiskCache = NULL, chunkSize = NULL,
  cacheDir = ".",
  cacheBase = ".blockConsModsCache",

  # Behaviour
  collectGarbage = FALSE,
  verbose = 1, indent = 0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individualData}{
Individual data from which the consensus is to be calculated. It can be either a list or a
\code{\link{multiData}} structure. Each element in \code{individulData} can in turn either be a numeric
obeject (vector, matrix or array) or a \code{\link{BlockwiseData}} structure.
}
  \item{consensusOptions}{
A list of class \code{ConsensusOptions} that contains options for the consensus calculation. A suitable list
can be obtained by calling function \code{\link{newConsensusOptions}}.
}
  \item{useBlocks}{
When \code{individualData} contains \code{\link{BlockwiseData}}, this argument can be an
integer vector with indices of blocks for which the calculation should be performed.
}
  \item{randomSeed}{
If non-\code{NULL}, the function will save the current state of the random generator, set the given seed,
and restore the random seed to its original state upon exit. If \code{NULL}, the seed is not set nor is it
restored on exit.
}
  \item{saveCalibratedIndividualData}{
Logical: should calibrated individual data be saved?
}
  \item{calibratedIndividualDataFilePattern}{
Pattern from which file names for saving calibrated individual data are determined. The conversions
\code{\%a}, \code{\%s} and \code{\%b} will be replaced by analysis name, set number and block number,
respectively.}
  \item{saveConsensusData}{
Logical: should final consensus be saved (\code{TRUE}) or returned in the return value (\code{FALSE})?
If \code{NULL}, data will be saved only if input data were blockwise data saved on disk rather than held in memory
}
  \item{consensusDataFileNames}{
Pattern from which file names for saving the final consensus are determined. The conversions
\code{\%a} and \code{\%b} will be replaced by analysis name and block number,
respectively.}
  \item{getCalibrationSamples}{
When calibration method in the \code{consensusOptions} component of \code{ConsensusTree} is \code{"single
quantile"}, this logical argument determines whether the calibration samples should be retuned within the
return value.
}
  \item{useDiskCache}{
Logical: should disk cache be used for consensus calculations? The disk cache can be used to sture chunks of
calibrated data that are small enough to fit one chunk from each set into memory (blocks may be small enough
to fit one block of one set into memory, but not small enogh to fit one block from all sets in a consensus
calculation into memory at the same time). Using disk cache is slower but lessens the memry footprint of
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
  \item{collectGarbage}{
Logical: should garbage collection be forced after each major calculation?
}
  \item{verbose}{Integer level of verbosity of diagnostic messages. Zero means silent, higher values make the
output
progressively more and more verbose. }

  \item{indent}{Indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
Consensus is defined as the element-wise (also known as "parallel") quantile of the individual data at
probability given by the \code{consensusQuantile} element of \code{consensusOptions}. Depending on the value
of component \code{calibration} of \code{consensusOptions}, the individual data are first calibrated. For
\code{consensusOptions$calibration="full quantile"}, the individual data are quantile normalized using
\code{\link{normalize.quantiles}}. For
\code{consensusOptions$calibration="single quantile"}, the individual data are raised to a power such that
the quantiles at probability \code{consensusOptions$calibrationQuantile} are the same. 
For \code{consensusOptions$calibration="none"}, the individual data are not calibrated.
}
\value{
A list with the following components:
\item{consensusData}{A \code{\link{BlockwiseData}} list containing the consensus.}
\item{nSets}{Number of input data sets.}
\item{saveCalibratedIndividualData}{Copy of the input \code{saveCalibratedIndividualData}.}
\item{calibratedIndividualData}{If input \code{saveCalibratedIndividualData} is \code{TRUE}, 
   a list in which each component is a \code{\link{BlockwiseData}} structure containing the calibrated
individual data for the corresponding input individual data set.}
\item{calibrationSamples}{If \code{consensusOptions$calibration} is \code{"single quantile"} and 
\code{getCalibrationSamples} is \code{TRUE}, a list in which eahc component contains the calibration samples
for the corresponding input individual data set.}
\code{originCount}{A vector of length \code{nSets} that
contains, for each set, the number of (calibrated) elements that were less than or equal the consensus for that element.}
}
\references{
Consensus network analysis was originally described in 
Langfelder P, Horvath S. Eigengene networks for studying the relationships
between co-expression modules. BMC Systems Biology 2007, 1:54 
http://www.biomedcentral.com/1752-0509/1/54
}
\author{
Peter Langfelder
}
\seealso{
\code{\link[preprocessCore]{normalize.quantiles}} for quantile normalization.
}
\keyword{misc}
