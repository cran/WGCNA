\name{goodSamplesGenes}
\alias{goodSamplesGenes}
\title{ Iterative filtering of samples and genes with too many missing entries }
\description{
   This function checks data for missing entries and zero-variance genes, 
and returns a list of samples and genes 
that pass criteria maximum number of missing values. If necessary, the filtering is iterated.

}
\usage{
goodSamplesGenes(
  datExpr, 
  minFraction = 1/2, 
  minNSamples = ..minNSamples, 
  minNGenes = ..minNGenes, 
  tol = NULL,
  verbose = 1, indent = 0)
}
\arguments{
  \item{datExpr}{  expression data. A data frame in which columns are genes and rows ar samples. }

  \item{minFraction}{ minimum fraction of non-missing samples for a gene to be considered good. }

  \item{minNSamples}{ minimum number of non-missing samples for a gene to be considered good.  }

  \item{minNGenes}{ minimum number of good genes for the data set to be considered fit for analysis. If
the actual number of good genes falls below this threshold, an error will be issued. }

  \item{tol}{ an optional 'small' number to compare the variance against. Defaults to the square of
 \code{1e-10 * max(abs(datExpr), na.rm = TRUE)}. The reason of comparing the variance to this number, rather
than
zero, is that the fast way of computing variance used by this function sometimes causes small numerical 
overflow errors which make variance of constant vectors slightly non-zero; comparing the variance to
\code{tol} rather than zero prevents the retaining of such genes as 'good genes'.}

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
  This function iteratively identifies samples and genes with too many missing entries and genes with
zero variance. Iterations may be
required since excluding samples effectively changes criteria on genes and vice versa. The process is
repeated until the lists of good samples and genes are stable.
The constants \code{..minNSamples} and \code{..minNGenes} are both set to the value 4.

}
\value{

  A list with the foolowing components:
  \item{goodSamples}{   A logical vector with one entry per sample that is \code{TRUE} if the sample is
considered good and \code{FALSE} otherwise.  }

  \item{goodGenes}{  A logical vector with one entry per gene that is \code{TRUE} if the gene is
considered good and \code{FALSE} otherwise.  }

}
\author{ Peter Langfelder }
\seealso{ \code{\link{goodSamples}}, \code{\link{goodGenes}} }

\keyword{ misc }% __ONLY ONE__ keyword per line
