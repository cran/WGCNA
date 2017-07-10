\name{bicovWeights}
\alias{bicovWeights}
\alias{bicovWeightFactors}
\alias{bicovWeightsFromFactors}
\title{
Weights used in biweight midcovariance
}
\description{
Calculation of weights and the intermediate weight factors 
used in the calculation of biweight midcovariance and midcorrelation. The
weights are designed such that outliers get smaller weights; the weights become zero for data points more than
9 median absolute deviations from the median.
}
\usage{
bicovWeights(
   x, 
   pearsonFallback = TRUE, 
   maxPOutliers = 1,
   outlierReferenceWeight = 0.5625,
   defaultWeight = 0)

bicovWeightFactors(
   x, 
   pearsonFallback = TRUE, 
   maxPOutliers = 1, 
   outlierReferenceWeight = 0.5625,
   defaultFactor = NA)

bicovWeightsFromFactors(
   u, 
   defaultWeight = 0)
}
\arguments{
  \item{x}{
A vector or a two-dimensional array (matrix or data frame). 
If two-dimensional, the weights will be calculated separately on each column.
}

  \item{u}{
A vector or matrix of weight factors, usually calculated by \code{bicovWeightFactors}.
}
  \item{pearsonFallback}{
Logical: if the median absolute deviation is zero, should standard deviation be substituted?
}
  \item{maxPOutliers}{
Optional specification of the maximum proportion of outliers, i.e., data with weights equal to
\code{outlierReferenceWeight} below.
}
  \item{outlierReferenceWeight}{A number between 0 and 1 specifying what is to be considered an outlier when
calculating the proportion of outliers.}

  \item{defaultWeight}{Value used for weights that correspond to a finite \code{x} but the weights themselves
would not be finite, for example, when a column in \code{x} is constant.}

  \item{defaultFactor}{Value used for factors that correspond to a finite \code{x} but the weights themselves
would not be finite, for example, when a column in \code{x} is constant.}

}
\details{
These functions are based on Equations (1) and (3) in Langfelder and Horvath (2012). The weight factor is denoted
\code{u} in that article.

Langfelder and Horvath (2012) also describe the Pearson fallback and 
maximum proportion of outliers in detail. For a full
discussion of the biweight midcovariance and midcorrelation, see Wilcox (2005).
}

\value{
  A vector or matrix of the same dimensions as the input \code{x} giving the bisquare weights
(\code{bicovWeights} and \code{bicovWeightsFromFactors}) or the bisquare factors
(\code{bicovWeightFactors}).
}
\references{
Langfelder P, Horvath S (2012) Fast R Functions for Robust Correlations and Hierarchical Clustering Journal of
Statistical Software 46(11) 1-17 PMID: 23050260 PMCID: PMC3465711 
Wilcox RR (2005). Introduction to Robust Estimation and Hypothesis Testing. 2nd edition.
Academic Press, Section 9.3.8, page 399 as well as Section 3.12.1, page 83. 
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{bicor}}
}
\examples{
x = rnorm(100);
x[1] = 10;
plot(x, bicovWeights(x));
}

\keyword{misc}
