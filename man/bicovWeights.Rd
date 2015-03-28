\name{bicovWeights}
\alias{bicovWeights}
\title{
Weights used in biweight midcovariance
}
\description{
The function calculates weights used in the calculation of biweight midcovariance and midcorrelation. The
weights are designed such that outliers get smaller weights; the weights become zero for data points more than
9 median absolute deviations from the median.
}
\usage{
bicovWeights(x, pearsonFallback = TRUE, maxPOutliers = 1)
}
\arguments{
  \item{x}{
A vector or a two-dimensional array (matrix or data frame). 
If two-dimensional, the weights will be calculated separately on each column.
}
  \item{pearsonFallback}{
Logical: if the median absolute deviation is zero, should standard deviation be substituted?
}
  \item{maxPOutliers}{
Optional specification of the maximum proportion of outliers, i.e., data with weights equal 0. 
}
}
\value{
  A vector or matrix of the same dimensions as the input \code{x} giving the weights.
}
\references{
This function is based on Equation (3) in

Langfelder P, Horvath S (2012) Fast R Functions for Robust Correlations and Hierarchical Clustering Journal of
Statistical Software 46(11) 1-17 PMID: 23050260 PMCID: PMC3465711 

That article also describes the Pearson fallback and maximum proportion of outliers in detail. For a full
discussion of the biweight midcovariance and midcorrelation, see

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
