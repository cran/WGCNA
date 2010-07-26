\name{rankPvalue}
\alias{rankPvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the p-value for ranking consistently high (or low) on multiple lists
}
\description{
The function rankPvalue calculates the p-value for observing that an object (corresponding to a row of the input
data frame \code{datS}) has a consistently high ranking (or low ranking) according to multiple ordinal scores
(corresponding to the columns of the input data frame \code{datS}).
}
\usage{
rankPvalue(datS, columnweights = NULL, 
           na.last = "keep", ties.method = "average", 
           calculateQvalue = TRUE, pValueMethod = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datS}{
a data frame whose rows represent objects that will be ranked. Each column of \code{datS} represents an
ordinal variable (which can take on negative values). The columns correspond to (possibly signed) object
significance measures, e.g., statistics (such as Z statistics), ranks, or correlations. 
}
  \item{columnweights}{
allows the user to input a vector of non-negative numbers reflecting weights for the different columns of
\code{datZ}.
If it is set to \code{NULL} then all weights are equal.
}
  \item{na.last}{
controls the treatment of missing values (NAs) in the rank function. If \code{TRUE}, missing values in the data are
put last (i.e. they get the highest rank values). If \code{FALSE}, they are put first; 
if \code{NA}, they are removed; if
\code{"keep"} they are kept with rank NA. See \code{\link{rank}} for more details.
}
  \item{ties.method}{
represents the ties method used in the rank function for the percentile rank method. See \code{\link{rank}} for
more details.}

  \item{calculateQvalue}{ logical: should q-values be calculated? If set to TRUE then the function calculates
corresponding q-values (local false discovery rates) using the qvalue package, see Storey JD and Tibshirani R.
(2003). This option assumes that qvalue package has been installed. }

  \item{pValueMethod}{
determines which method is used for calculating p-values. By default it is set to "all", i.e. both methods are
used. If it is set to "rank" then only the percentile rank method is used. If it set to "scale" then only the
scale method will be used.
}
}
\details{
The function calculates asymptotic p-values (and optionally q-values) for testing the null hypothesis that the
values in the columns of datS are independent. This allows us to find objects (rows) with consistently high (or
low) values across the columns.

Example: Imagine you have 5 vectors of Z statistics corresponding to the columns of datS. Further assume that a
gene has ranks 1,1,1,1,20 in the 5 lists. It seems very significant that the gene ranks number 1 in 4 out of the
5 lists. The function rankPvalue can be used to calculate a p-value for this occurrence.

The function uses the central limit theorem to calculate asymptotic p-values for two types of test statistics
that measure consistently high or low ordinal values.
The first method (referred to as percentile rank method) leads to accurate estimates of p-values if datS has at
least 4 columns but it can be overly conservative. 
The percentile rank method replaces each column datS by the ranked version rank(datS[,i]) (referred to ask low
ranking) and by rank(-datS[,i]) (referred to as high ranking). Low ranking and high ranking allow one to find
consistently small values or  consistently large values of datS, respectively.  All ranks are divided by the
maximum rank so that the result lies in the unit interval [0,1]. In the following, we refer to rank/max(rank) as
percentile rank. For a given object (corresponding to a row of datS) the observed percentile rank follows
approximately a uniform distribution under the null hypothesis. The test statistic is defined as the sum of the
percentile ranks (across the columns of datS). Under the null hypothesis that there is no relationship between
the rankings of the columns of datS, this (row sum) test statistic follows a distribution that is given by the
convolution of random uniform distributions. Under the null hypothesis, the individual percentile ranks are
independent and one can invoke the central limit theorem to argue that the row sum test statistic follows
asymptotically a normal distribution.  It is well-known that the speed of convergence to the normal distribution
is extremely fast in case of  identically distributed uniform distributions. Even when datS has only  4 columns,
the difference between the normal approximation and the exact distribution is negligible in practice (Killmann
et al 2001).
In summary, we use the central limit theorem to argue that the sum of the percentile ranks follows a normal
distribution whose mean and variance can be calculated using the fact that the mean value of a uniform random
variable (on the unit interval) equals 0.5 and its variance equals 1/12.

The second method for calculating p-values is referred to as scale method. It is often more powerful but its
asymptotic p-value can only be trusted if either datS has a lot of columns or if the ordinal scores (columns of
datS)  follow an approximate normal distribution.  The scale method scales (or standardizes) each ordinal
variable (column of datS) so that it has mean 0 and variance 1. Under the null hypothesis of independence, the
row sum follows approximately a normal distribution if the assumptions of  the central limit theorem are met.
In practice, we find that the second approach is often more powerful but it makes more distributional
assumptions (if datS has few columns).
}
\value{

A list whose actual content depends on which p-value methods is selected, and whether q0values are calculated.
The following inner components are calculated, organized in outer components \code{datoutrank} and
\code{datoutscale},:

\item{pValueExtremeRank}{This is the minimum between pValueLowRank and
pValueHighRank, i.e. min(pValueLow, pValueHigh)}

\item{pValueLowRank}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the rank method.}

\item{pValueHighRank}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the rank method.}

\item{pValueExtremeScale}{This is the minimum between pValueLowScale and
pValueHighScale, i.e. min(pValueLow, pValueHigh)}

\item{pValueLowScale}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the Scale method.}

\item{pValueHighScale}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the Scale method.}

\item{qValueExtremeRank}{local false discovery rate (q-value) corresponding to the p-value
pValueExtremeRank}

\item{qValueLowRank}{local false discovery rate (q-value) corresponding to the p-value
pValueLowRank}

\item{qValueHighRank}{local false discovery rate (q-value) corresponding to the p-value
pValueHighRank}

\item{qValueExtremeScale}{local false discovery rate (q-value) corresponding to the p-value
pValueExtremeScale}

\item{qValueLowScale}{local false discovery rate (q-value) corresponding to the p-value
pValueLowScale}

\item{qValueHighScale}{local false discovery rate (q-value) corresponding to the p-value
pValueHighScale}

}
\references{
Killmann F,  VonCollani E (2001) A Note on the Convolution of the Uniform and Related Distributions and Their
Use in Quality Control. Economic Quality Control Vol 16 (2001), No. 1, 17-41.ISSN 0940-5151

Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide experiments. Proceedings of the
National Academy of Sciences, 100: 9440-9445. }
\author{
Steve Horvath
}

\seealso{
\code{\link{rank}}, \code{\link{qvalue}}
}
\keyword{ misc}% __ONLY ONE__ keyword per line
