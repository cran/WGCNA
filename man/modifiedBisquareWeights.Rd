\name{modifiedBisquareWeights}
\alias{modifiedBisquareWeights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modified Bisquare Weights 
}
\description{
Calculation of bisquare weights and the intermediate weight factors
similar to those used in the calculation of biweight midcovariance and midcorrelation. The
weights are designed such that outliers get smaller weights; the weights become zero for data points more than
9 median absolute deviations from the median.
}

\usage{
modifiedBisquareWeights(
  x,
  removedCovariates = NULL,
  pearsonFallback = TRUE,
  maxPOutliers = 0.05,
  outlierReferenceWeight = 0.1,
  groupsForMinWeightRestriction = NULL,
  minWeightInGroups = 0,
  maxPropUnderMinWeight = 1,
  defaultWeight = 1,
  getFactors = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix of numeric observations with variables (features) in columns and observations (samples) in rows. 
Weights will be calculated separately for each column.
}

  \item{removedCovariates}{
Optional matrix or data frame of variables that are to be regressed out of each column
of \code{x} before calculating the weights. If given, must have the same number of rows as \code{x}. }

  \item{pearsonFallback}{
Logical: for columns of \code{x} that have zero median absolute deviation (MAD), should the appropriately scaled standard
deviation be used instead?}

  \item{maxPOutliers}{
Optional numeric scalar between 0 and 1. Specifies the maximum proportion of outliers in each column, 
i.e., data with weights equal to
\code{outlierReferenceWeight} below.  }

  \item{outlierReferenceWeight}{A number between 0 and 1 specifying what is to be considered an outlier when
calculating the proportion of outliers.}

  \item{groupsForMinWeightRestriction}{
An optional vector with length equal to the number of samples (rows) in \code{x} giving a categorical variable. The
output factors and weights are adjusted such that in samples at each level of the variable, the weight is below
\code{minWeightInGroups} in a fraction of samples that is at most \code{maxPropUnderMinWeight}.}

  \item{minWeightInGroups}{
A threshold weight, see \code{groupsForMinWeightRestriction} and details.
}
  \item{maxPropUnderMinWeight}{
A proportion (number between 0 and 1). See  \code{groupsForMinWeightRestriction} and details.
}
  \item{defaultWeight}{Value used for weights that would be undefined or not finite, for example, when a 
        column in \code{x} is constant.}

  \item{getFactors}{
Logical: should the intermediate weight factors be returned as well?  }
}

\details{
Weights are calculated independently for each column of \code{x}. Denoting a column of \code{x} as \code{y}, the weights
are calculated as \eqn{(1-u^2)^2}{(1-u^2)^2} where \code{u} is defined as 
\eqn{\min(1, |y-m|/(9MMAD))}{min(1, abs(y-m)/(9 * MMAD))}. Here \code{m} is the median
of the column \code{y} and \code{MMAD} is the modified median absolute deviation. We call the expression 
\eqn{|y-m|/(9 MMAD)}{abs(y-m)/(9 * MMAD)} 
the weight factors. Note that outliers are observations with high (>1) weight factors for outliers but low (zero) weights.

The calculation of \code{MMAD} starts
with calculating the (unscaled) median absolute deviation of the column \code{x}. If the median absolute deviation is
zero and \code{pearsonFallback} is TRUE, it is replaced by the standard deviation 
(multiplied by \code{qnorm(0.75)} to make it asymptotically consistent with
MAD). The following two conditions are then
checked: (1) The proportion of weights below \code{outlierReferenceWeight} is at most \code{maxPOutliers} 
and (2) if \code{groupsForMinWeightRestriction} has non-zero length, then for each individual level in
\code{groupsForMinWeightRestriction}, the proportion of samples with weights less than \code{minWeightInGroups} is at
most \code{maxPropUnderMinWeight}. (If \code{groupsForMinWeightRestriction} is zero-length, the second condition is
considered trivially satisfied.) If both conditions are met, \code{MMAD} equals the median absolute deviation, MAD. If
either condition is not met, MMAD equals the lowest number for which both conditions are met.   
}
\value{
When the input \code{getFactors} is \code{TRUE}, a list with two components:
\item{weights}{A matrix of the same dimensions and \code{dimnames} as the input \code{x} giving the weights of the
individual observations in \code{x}.}
\item{factors}{A matrix of the same form as \code{weights} giving the weight factors.}
When the input \code{getFactors} is \code{FALSE}, the function returns the matrix of weights.
}
\references{

A full description of the weight calculation can be found, e.g., in Methods section of 

Wang N, Langfelder P, et al (2022) 
Mapping brain gene coexpression in daytime transcriptomes unveils diurnal 
molecular networks and deciphers perturbation gene signatures. Neuron. 2022 Oct 19;110(20):3318-3338.e9. 
PMID: 36265442; PMCID: PMC9665885.
\doi{10.1016/j.neuron.2022.09.028}

Other references include, in reverse chronological order, 

Peter Langfelder, Steve Horvath (2012)
Fast R Functions for Robust Correlations and Hierarchical Clustering.
Journal of Statistical Software, 46(11), 1-17.
\url{https://www.jstatsoft.org/v46/i11/}

"Introduction to Robust Estimation and Hypothesis Testing", Rand Wilcox, Academic Press, 1997.

"Data Analysis and Regression: A Second Course in Statistics", Mosteller and Tukey, Addison-Wesley,
 1977, pp. 203-209.

}
\author{
Peter Langfelder
}

\seealso{
\code{bicovWeights} for a simpler, less flexible calculation. 
}
\keyword{misc}
