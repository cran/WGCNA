\name{hierarchicalConsensusCalculation}
\alias{hierarchicalConsensusCalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical consensus calculation
}
\description{
Hierarchical consensus calculation with optional data calibration.
}
\usage{
hierarchicalConsensusCalculation(
   individualData,

   consensusTree,

   level = 1,
   useBlocks = NULL,
   randomSeed = NULL,
   saveCalibratedIndividualData = FALSE,
   calibratedIndividualDataFilePattern = 
         "calibratedIndividualData-\%a-Set\%s-Block\%b.RData",

   # Return options: the data can be either saved or returned but not both.
   saveConsensusData = TRUE,
   consensusDataFileNames = "consensusData-\%a-Block\%b.RData",
   getCalibrationSamples= FALSE,

   # Return the intermediate results as well?
   keepIntermediateResults = FALSE,

   # Internal handling of data
   useDiskCache = NULL, 
   chunkSize = NULL,
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
object (vector, matrix or array) or a \code{\link{BlockwiseData}} structure.
}
  \item{consensusTree}{
A list specifying the consensus calculation. See details.
}
  \item{level}{
Integer which the user should leave at 1.  This serves to keep default set names unique.
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
}
  \item{consensusDataFileNames}{
Pattern from which file names for saving the final consensus are determined. The conversions
\code{\%a} and \code{\%b} will be replaced by analysis name and block number,
respectively.}
  \item{getCalibrationSamples}{
When calibration method in the \code{consensusOptions} component of \code{ConsensusTree} is \code{"single
quantile"}, this logical argument determines whether the calibration samples should be returned within the
return value.
}
  \item{keepIntermediateResults}{
Logical: should results of intermediate consensus calculations (if any) be kept? These are always returned
as \code{BlockwiseData} whose data are saved to disk.
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
  \item{collectGarbage}{
Logical: should garbage collection be forced after each major calculation?
}
  \item{verbose}{Integer level of verbosity of diagnostic messages. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{Indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
This function calculates consensus in a hierarchical manner, using a separate (and possibly different) set of
consensus options at each step. The "recipe" for the consensus calculation is supplied in the argument
\code{consensusTree}. 

The argument \code{consensusTree} should have the following components: (1) \code{inputs} must be either a
character vector whose components match \code{names(inputData)}, or consensus trees in the own right.
(2) \code{consensusOptions} must be a list of class \code{"ConsensusOptions"} that specifies options for
calculating the consensus. A suitable set of options can be obtained by calling
\code{\link{newConsensusOptions}}. (3) Optionally, the component \code{analysisName} can be a single
character string giving the name for the analysis. When intermediate results are returned, they are returned
in a list whose names will be set from \code{analysisName} components, if they exist.

The actual consensus calculation at each level of the consensus tree
is carried out in function \code{\link{consensusCalculation}}. The consensus options for each individual
consensus calculation are independent from one another, i.e., the consensus options for different steps can
be different.
}
\value{
A list containing the output of the top level call to \code{\link{consensusCalculation}}; if
\code{keepIntermediateResults} is \code{TRUE}, component \code{inputs} contains a (possibly recursive) list
of the results of intermediate consensus calculations. Names of the \code{inputs} list are taken from the
corresponding \code{analysisName} components if they exist, otherwise from names of the corresponding
\code{inputs} components of the supplied \code{consensusTree}.  See example below for an example of a
relatively simple consensus tree.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{newConsensusOptions}} for obtaining a suitable list of consensus options;

\code{\link{consensusCalculation}} for the actual calculation of a consensus that underpins this function.
}
\examples{
# We generate 3 simple matrices
set.seed(5)
data = replicate(3, matrix(rnorm(10*100), 10, 100))
names(data) = c("Set1", "Set2", "Set3");
# Put together a consensus tree. In this example the final consensus uses 
# as input set 1 and a consensus of sets 2 and 3. 

# First define the consensus of sets 2 and 3:
consTree.23 = newConsensusTree(
           inputs = c("Set2", "Set3"),
           consensusOptions = newConsensusOptions(calibration = "none",
                               consensusQuantile = 0.25),
           analysisName = "Consensus of sets 1 and 2");

# Now define the final consensus
consTree.final = newConsensusTree(
   inputs = list("Set1", consTree.23),
   consensusOptions = newConsensusOptions(calibration = "full quantile",
                               consensusQuantile = 0),
   analysisName = "Final consensus");

consensus = hierarchicalConsensusCalculation(
  individualData = data,
  consensusTree = consTree.final,
  saveConsensusData = FALSE,
  keepIntermediateResults = FALSE)

names(consensus)
}
\keyword{misc}
