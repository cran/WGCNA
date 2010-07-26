\name{goodGenes}
\alias{goodGenes}
\title{ Filter genes with too many missing entries }
\description{
  This function checks data for missing entries and returns a list of genes that have non-zero variance
and pass two criteria on maximum
number of missing values: the fraction of missing values must be below a given threshold and the total number
of missing samples must be below a given threshold. 
}
\usage{
goodGenes(datExpr, 
          useSamples = NULL, 
          useGenes = NULL, 
          minFraction = 1/2, 
          minNSamples = ..minNSamples, 
          minNGenes = ..minNGenes, 
          verbose = 1, indent = 0)
}
\arguments{
  \item{datExpr}{  expression data. A data frame in which columns are genes and rows ar samples. }

  \item{useSamples}{ optional specifications of which samples to use for the check. Should be a logical
vector; samples whose entries are \code{FALSE} will be ignored for the missing value counts. Defaults to
using all samples.}

  \item{useGenes}{ optional specifications of genes for which to perform the check. Should be a logical 
vector; genes whose entries are \code{FALSE} will be ignored. Defaults to 
using all genes.}

  \item{minFraction}{ minimum fraction of non-missing samples for a gene to be considered good. }

  \item{minNSamples}{ minimum number of non-missing samples for a gene to be considered good.  }

  \item{minNGenes}{ minimum number of good genes for the data set to be considered fit for analysis. If
the actual number of good genes falls below this threshold, an error will be issued. }

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
  A logical vector with one entry per gene that is \code{TRUE} if the gene is considered good and
\code{FALSE} otherwise. Note that all genes excluded by \code{useGenes} are automatically assigned
\code{FALSE}. 
}
\author{ Peter Langfelder and Steve Horvath }
\seealso{ \code{\link{goodSamples}}, \code{\link{goodSamplesGenes}} }
\keyword{ misc }% __ONLY ONE__ keyword per line
