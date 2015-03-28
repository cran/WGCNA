\name{empiricalBayesLM}
\alias{empiricalBayesLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Empirical Bayes-moderated adjustment for unwanted covariates
}
\description{
This functions removes variation in high-dimensional data 
due to unwanted covariates while preserving variation due to retained covariates. To prevent numerical
instability, it uses Empirical bayes-moderated linear regression, optionally in a robust (outlier-resistant)
form.
}
\usage{
empiricalBayesLM(
  data,
  retainedCovariates,
  removedCovariates,
  weights = NULL,
  weightType = c("apriori", "empirical"),
  stopOnSmallWeights = TRUE,
  tol = 1e-4, maxIterations = 1000,
  scaleMeanToSamples = NULL,
  robustPriors = FALSE,
  automaticWeights = c("none", "bicov"),
  aw.maxPOutliers = 0.1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A 2-dimensional matrix or data frame of numeric data to be adjusted. Variables (for example, genes or
methylation profiles) should be in columns and observations (samples) should be in rows.
}
  \item{retainedCovariates}{
A vector or two-dimensional object (matrix or data frame) giving the covariates whose effect on the data is to be retained. May be
\code{NULL} if there are no such "retained" covariates.
}
  \item{removedCovariates}{
A vector or two-dimensional object (matrix or data frame) giving the covariates whose effect on the data is to
be removed. At least one such covariate must be given.
}
  \item{weights}{
Optional 2-dimensional matrix or data frame of the same dimensions as \code{data} giving weights for each
entry in \code{data}.
}
  \item{weightType}{
One of (unique abbreviations of) \code{"apriori"} or \code{"empirical"}. Determines whether a standard
(\code{"apriori"}) or a modified (\code{"empirical"}) weighted regression is used. The \code{"apriori"} choice is
suitable for weights that have been determined without knowledge of the actual \code{data}, while
\code{"empirical"} is appropriate for situations where one wants to down-weigh cartain entries of \code{data}
because they may be outliers. In either case, the weights should be determined in a way that is independent of 
the covariates (both retained and removed).
}
  \item{stopOnSmallWeights}{
Logical: should presence of small \code{"apriori"} weights trigger an error? Because standard weighted regression
assumes that all weights are non-zero (otherwise estimates of standard errors will be biased), this function
will by default complain about the presence of too small \code{"apriori"} weights.
}
  \item{tol}{
Convergence criterion used in the numerical equation solver. When the relative change in coefficients falls
below this threshold, the system will be considered to have converged.
}
  \item{maxIterations}{
Maximum number of iterations to use.
}
  \item{scaleMeanToSamples}{
Optional specification of samples (given as a vector of indices) to whose means the resulting adjusted data
should be scaled (more precisely, shifted). If not given, the mean of all samples will be used.
}
  \item{robustPriors}{
Logical: should robust priors be used? This essentially means replacing mean by median and covariance by
biweight mid-covariance.
}
  \item{automaticWeights}{
One of (unique abrreviations of) \code{"none"} or \code{"bicov"}, instructing the function to calculate
weights from the given \code{data}. Value \code{"none"} will result in trivial weights; value \code{"bicov"}
will result in biweight midcovariance weights being used.
}
  \item{aw.maxPOutliers}{
If \code{automaticWeights} above is \code{"bicov"}, this argument gets passed to the function
\code{\link{bicovWeights}} and determines the maximum proportion of outliers in calculating the weights. See
\code{\link{bicovWeights}} for more details.
}
}
\details{
This function uses Empirical Bayes-moderated (EB) linear regression to remove variation in \code{data} due to the
variables in \code{removedCovariates} while retaining variation due to variables in \code{retainedCovariates},
if any are given. The EB step uses simple normal priors on the regression coefficients and inverse gamma 
priors on the
variances. The procedure starts with multivariate ordinary linear regression of individual columns in
\code{data} on \code{retainedCovariates} and \code{removedCovariates}. To make the coefficients comparable,
columns of \code{data} are scaled to (weighted if weights are given) mean 0 and variance 1.
The resulting regression coefficients are used to
determine the parameters of the normal prior (mean, covariance, and inverse gamma or median and biweight
mid-covariance if robust priors are used), and the variances are used to determine the parameters of the
inverse gamma prior. The EB step then essentially shrinks the coefficients toward their means, with the amount
of shrinkage determined by the prior covariance.

Using appropriate weights can make the data adjustment robust to outliers. This can be achieved automatically
by using the argument \code{automaticWeights = "bicov"}. When bicov weights are used, we also recommend
setting the argument \code{maxPOutliers} to a maximum proportion of samples that could be outliers. This is
especially important if some of the design variables are binary and can be expected to have a strong effect on
some of the columns in \code{data}, since standard biweight midcorrelation (and its weights) do not work well
on bimodal data. 

The automatic bicov weights are determined from \code{data} only. It is implicitly assumed that there are no
outliers in the retained and removed covariates. Outliers in the covariates are more difficult to work with
since, even if the regression is made robust to them, they can influence the adjusted values for the sample in
which they appear. Unless the the covariate outliers can be attributed to a relevant vriation in experimental
conditions, samples with covariate outliers are best removed entirely before calling this function.

}
\value{
A list with the following components:
\item{adjustedData}{A matrix of the same dimensions as the input \code{data}, giving the adjusted data. If
input \code{data} has non-NULL \code{dimnames}, these are copied.}

\item{residuals}{A matrix of the same dimensions as the input \code{data}, giving the residuals,
that is, adjusted data with zero means.}

\item{coefficients}{A matrix of regression coefficients. Rows correspond to the design matrix variables 
(mean, retained and removed covariates) and columns correspond to the variables (columns) in \code{data}.}

\item{coefficiens.scaled}{A matrix of regression coefficients corresponding to columns in \code{data} scaled
to mean 0 and variance 1.}

\item{sigmaSq}{Estimated error variances (one for each column of input \code{data}.}

\item{sigmaSq.scaled}{Estimated error variances corresponding to columns in \code{data} scaled
to mean 0 and variance 1.}

\item{fittedValues}{Fitted values calculated from the means and coefficients corresponding to the removed
covariates, i.e., roughly the values that are subtracted out of the data.}

\item{adjustedData.OLS}{A matrix of the same dimensions as the input \code{data}, giving the data adjusted by
ordinary least squares. This component should only be used for diagnostic purposes, not as input for further
downstream analyses, as the OLS adjustment is inferior to EB adjustment. }

\item{residuals.OLS}{A matrix of the same dimensions as the input \code{data}, giving the residuals obtained
from ordinary least squares regression, that is, OLS-adjusted data with zero means.}

\item{coefficients.OLS}{A matrix of ordinary least squares regression coefficients. 
Rows correspond to the design matrix variables        
(mean, retained and removed covariates) and columns correspond to the variables (columns) in \code{data}.}

\item{coefficiens.OLS.scaled}{A matrix of ordinary least squares regression coefficients corresponding to columns
in \code{data} scaled to mean 0 and variance 1.  These coefficients are used to calculate priors for the EB step.}

\item{sigmaSq.OLS}{Estimated OLS error variances (one for each column of input \code{data}.}

\item{sigmaSq.OLS.scaled}{Estimated OLS error variances corresponding to columns in \code{data} scaled
to mean 0 and variance 1. These are used to calculate variance priors for the EB step.}

\item{fittedValues.OLS}{OLS fitted values calculated from the means and coefficients corresponding to the removed
covariates.}


}
\author{
Peter Langfelder
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bicovWeights}} for suitable weights that make the adjustment robust to outliers.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
