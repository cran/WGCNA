\name{newConsensusOptions}
\alias{newConsensusOptions}
\alias{ConsensusOptions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a list holding consensus calculation options.
}
\description{
This function creates a list of class \code{ConsensusOptions} that holds options for consensus calculations.
This list holds options for a single-level analysis.
}
\usage{
newConsensusOptions(
      calibration = c("full quantile", "single quantile", "none"),

      # Simple quantile scaling options
      calibrationQuantile = 0.95,
      sampleForCalibration = TRUE, 
      sampleForCalibrationFactor = 1000,

      # Consensus definition
      consensusQuantile = 0,
      useMean = FALSE,
      setWeights = NULL,
      # Name to prevent clashing of files
      analysisName = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calibration}{
Calibration method. One of \code{"full quantile", "single quantile", "none"}
(or a unique abbreviation of one of them).}
  \item{calibrationQuantile}{
if \code{calibration} is \code{"single quantile"},
input data to a consensus calculation 
 will be scaled such that their \code{calibrationQuantile} quantiles will agree. }

  \item{sampleForCalibration}{ if \code{TRUE}, calibration quantiles will be determined from a sample of
network
similarities. Note that using all data can double the memory footprint of the function and the function
may fail. }

  \item{sampleForCalibrationFactor}{ Determines the number of samples for calibration: the number is
 \code{1/calibrationQuantile * sampleForCalibrationFactor}. Should be set well above 1 to ensure accuracy of
the sampled quantile. }

  \item{consensusQuantile}{Quantile at which consensus is to be defined. See details. 
 }
  \item{useMean}{
Logical: should the consensus be calculated using (weighted) mean rather than a quantile?
}
  \item{setWeights}{
Optional specification of weights when \code{useMean} is \code{TRUE}.
}
  \item{analysisName}{
Optional character string naming the consensus analysis. Useful for identifying partial consensus calculation
in hierarchical consensus analysis.
}
}

\value{
A list of type \code{ConsensusOptions} that holds copies of the input arguments.
}
\author{
Peter Langfelder
}
\keyword{misc}% __ONLY ONE__ keyword per line
