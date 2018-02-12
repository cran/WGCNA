\name{goodGenesMS}
\alias{goodGenesMS}
\title{Filter genes with too many missing entries across multiple sets }
\description{
  This function checks data for missing entries and returns a list of genes that have non-zero variance
in all sets and pass two criteria on
maximum number of missing values in each given set: the fraction of missing values must be below a given
threshold and the total number of missing samples must be below a given threshold. If weights are given,
entries whose relative weight is below a threshold will be considered missing. 

}
\usage{
goodGenesMS(
  multiExpr,
  multiWeights = NULL,
  useSamples = NULL,
  useGenes = NULL,
  minFraction = 1/2,
  minNSamples = ..minNSamples,
  minNGenes = ..minNGenes,
  tol = NULL,
  minRelativeWeight = 0.1,
  verbose = 1, indent = 0)

}
\arguments{

  \item{multiExpr}{  expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{multiWeights}{ optional observation weights in the same format (and dimensions) as \code{multiExpr}.}

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

  \item{tol}{ an optional 'small' number to compare the variance against. For each set in \code{multiExpr}, 
the default value is \code{1e-10 * max(abs(multiExpr[[set]]$data), na.rm = TRUE)}. 
The reason of comparing the variance to this number, rather than
zero, is that the fast way of computing variance used by this function sometimes causes small numerical 
overflow errors which make variance of constant vectors slightly non-zero; comparing the variance to
\code{tol} rather than zero prevents the retaining of such genes as 'good genes'.}

  \item{minRelativeWeight}{ observations whose relative weight is below
this threshold will be considered missing. Here relative weight is weight divided by the maximum weight in
the column (gene).}

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }


}
\details{
  The constants \code{..minNSamples} and \code{..minNGenes} are both set to the value 4. 

If weights are given, entries whose relative
weight (i.e., weight divided by maximum weight in the column or gene) will be considered missing.

    For most data sets, the fraction of missing samples criterion will be much more stringent than the
absolute number of missing samples criterion.

}
\value{
  A logical vector with one entry per gene that is \code{TRUE} if the gene is considered good and
\code{FALSE} otherwise. Note that all genes excluded by \code{useGenes} are automatically assigned
\code{FALSE}.

}
\author{ Peter Langfelder }
\seealso{ 
\code{\link{goodGenes}}, \code{\link{goodSamples}}, \code{\link{goodSamplesGenes}} for cleaning
individual sets separately;

\code{\link{goodSamplesMS}},  \code{\link{goodSamplesGenesMS}} for additional cleaning of multiple data
sets together. }
 
\keyword{ misc }% __ONLY ONE__ keyword per line
