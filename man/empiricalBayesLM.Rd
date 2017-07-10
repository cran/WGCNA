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
  removedCovariates,
  retainedCovariates = NULL,

  initialFitFunction = NULL,
  initialFitOptions = NULL,
  initialFitRequiresFormula = NULL,
  initialFit.returnWeightName = NULL,

  weights = NULL,
  weightType = c("apriori", "empirical"),
  stopOnSmallWeights = TRUE,
  tol = 1e-4, maxIterations = 1000,
  scaleMeanToSamples = NULL,
  robustPriors = FALSE,
  automaticWeights = c("none", "bicov"),
  aw.maxPOutliers = 0.1,
  minDesignDeviation = 1e-10,
  garbageCollectInterval = 50000,

  getOLSAdjustedData = TRUE,
  getResiduals = TRUE,
  getFittedValues = TRUE,
  getWeights = TRUE,
  getEBadjustedData = TRUE,

  verbose = 0, indent = 0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A 2-dimensional matrix or data frame of numeric data to be adjusted. Variables (for example, genes or
methylation profiles) should be in columns and observations (samples) should be in rows.
}
  \item{removedCovariates}{
A vector or two-dimensional object (matrix or data frame) giving the covariates whose effect on the data is to
be removed. At least one such covariate must be given.
}
  \item{retainedCovariates}{
A vector or two-dimensional object (matrix or data frame) giving the covariates whose effect on the data is
to be retained. May be \code{NULL} if there are no such "retained" covariates.
}

  \item{initialFitFunction}{
Function name to perform the initial fit. The default is to use the internal implementation of linear model
fitting. The function must take arguments \code{formula} and \code{data} or \code{x} and \code{y},
 plus possibly additional arguments.
The return value must be a list with component \code{coefficients}, either \code{scale} or
\code{residuals}, and weights must be returned in component specified by \code{initialFit.returnWeightName}. 
See \code{\link{lm}}, \code{\link{rlm}} and other standard fit functions for examples of
suitable functions.}

  \item{initialFitOptions}{
Optional specifications of extra arguments for \code{initialFitFunction}, apart from \code{formula} and
\code{data} or \code{x} and \code{y}. 
Defaults are provided for functions \code{\link[MASS]{rlm}} and \code{\link[robust]{lmRob}},
i.e., if either of these two functions are used as \code{initialFitFunction}, suitable initial fit options
will be chosen automatically.}

  \item{initialFitRequiresFormula}{
Logical: does the initial fit function need \code{formula} and \code{data} arguments? If \code{TRUE},
\code{initialFitFunction} will be called with arguments \code{formula} and \code{data}, otherwise with
arguments \code{x} and \code{y}.}

  \item{initialFit.returnWeightName}{
Name of the component of the return value of \code{initialFitFunction} that contains the weights used in the
fit. Suitable default value will be chosen automatically for \code{\link[MASS]{rlm}} and
\code{\link[robust]{lmRob}}.
}
  \item{weights}{
Optional 2-dimensional matrix or data frame of the same dimensions as \code{data} giving weights for each
entry in \code{data}. These weights will be used in the initial fit and are are separate from the ones returned by
\code{initialFitFunction} if it is specified.
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

  \item{minDesignDeviation}{
Minimum standard deviation for columns of the design matrix to be retained. Columns with standard deviations
below this number will be removed (effectively removing the corresponding terms from the design).
}

  \item{garbageCollectInterval}{
Number of variables after which to call garbage collection.
}

  \item{getOLSAdjustedData}{Logical: should data adjusted by ordinary least squares or by
\code{initialFitFunction}, if specified, be returned?} 

  \item{getResiduals}{Logical: should the residuals (adjusted values without the means) be returned?}

  \item{getFittedValues}{Logical: should fitted values be returned?}

  \item{getWeights}{Logical: should the final weights be returned?}

  \item{getEBadjustedData}{Logical: should the EB step be performed and the adjusted data returned? If this
is \code{FALSE}, the function acts as a rather slow but still potentially useful adjustment using standard
fit functions.}

  \item{verbose}{Level of verbosity. Zero means silent, higher values result in more diagnostic messages
being printed.}
  \item{indent}{Indentation of diagnostic messages. Each unit adds two spaces.}
}

\details{
This function uses Empirical Bayes-moderated (EB) linear regression to remove variation in \code{data} due to the
variables in \code{removedCovariates} while retaining variation due to variables in \code{retainedCovariates},
if any are given. The EB step uses simple normal priors on the regression coefficients and inverse gamma 
priors on the
variances. The procedure starts with multivariate ordinary linear regression of individual columns in
\code{data} on \code{retainedCovariates} and \code{removedCovariates}. Alternatively, the user may specify an
intial fit function (e.g., robust linear regression). To make the coefficients comparable,
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
A list with the following components (some of which may be missing depending on input options):
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

\item{weights}{A matrix of weights used in the regression models. The matrix has the same dimension as the
input \code{data}.}

\item{dataColumnValid}{Logical vector with one element per column of input \code{data}, indicating whether the
column was adjusted. Columns with zero variance or too many missing data cannot be adjusted.}

\item{dataColumnWithZeroVariance}{Logical vector with one element per column of input \code{data}, indicating
whether the  
column had zero variance.}

\item{coefficientValid}{Logical matrix of the dimension (number of covariates +1) times (number of
variables in \code{data}), indicating whether the corresponding regression coefficient is valid. Invalid
regression coefficients may be returned as missing values or as zeroes.}


}
\author{
Peter Langfelder
}

\seealso{
\code{\link{bicovWeights}} for suitable weights that make the adjustment robust to outliers.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
