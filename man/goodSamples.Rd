\name{goodSamples}
\alias{goodSamples}
\title{ Filter samples with too many missing entries }

\description{
 This function checks data for missing entries and returns a list of samples that pass two criteria on
maximum
number of missing values: the fraction of missing values must be below a given threshold and the total
number
of missing genes must be below a given threshold.

}
\usage{
goodSamples(
  datExpr, 
  weights = NULL,
  useSamples = NULL, 
  useGenes = NULL, 
  minFraction = 1/2, 
  minNSamples = ..minNSamples, 
  minNGenes = ..minNGenes, 
  minRelativeWeight = 0.1,
  verbose = 1, indent = 0)
}
\arguments{
  \item{datExpr}{  expression data. A data frame in which columns are genes and rows ar samples. }

  \item{weights}{optional observation weights in the same format (and dimensions) as \code{datExpr}.}

  \item{useSamples}{ optional specifications of which samples to use for the check. Should be a logical
vector; samples whose entries are \code{FALSE} will be ignored for the missing value counts. Defaults to
using all samples.}

  \item{useGenes}{ optional specifications of genes for which to perform the check. Should be a logical
vector; genes whose entries are \code{FALSE} will be ignored. Defaults to
using all genes.}

  \item{minFraction}{ minimum fraction of non-missing samples for a gene to be considered good. }

  \item{minNSamples}{ minimum number of good samples for the data set to be considered fit for analysis.
If the actual number of good samples falls below this threshold, an error will be issued. }

  \item{minNGenes}{ minimum number of non-missing samples for a sample to be considered good.  }

  \item{minRelativeWeight}{ observations whose weight divided by the maximum weight is below this threshold
will be considered missing. }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
  The constants \code{..minNSamples} and \code{..minNGenes} are both set to the value 4.
  For most data sets, the fraction of missing samples criterion will be much more stringent than the
absolute number of missing samples criterion.
}

\value{
    A logical vector with one entry per sample that is \code{TRUE} if the sample is considered good and
\code{FALSE} otherwise. Note that all samples excluded by \code{useSamples} are automatically assigned
\code{FALSE}.
}

\author{ Peter Langfelder and Steve Horvath }
\seealso{ \code{\link{goodSamples}}, \code{\link{goodSamplesGenes}} }
\keyword{ misc }% __ONLY ONE__ keyword per line
