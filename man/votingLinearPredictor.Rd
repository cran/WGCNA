\name{votingLinearPredictor}
\alias{votingLinearPredictor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Voting linear predictor }
\description{
Predictor based on univariate regression on all or selected given features 
that pools all predictions using weights
derived from the univariate linear models.
}
\usage{
votingLinearPredictor(
         x, y, xtest = NULL, 
         classify = FALSE, 
         CVfold = 0, 
         randomSeed = 12345, 
         assocFnc = "cor", assocOptions = "use = 'p'", 
         featureWeightPowers = NULL, priorWeights = NULL, 
         weighByPrediction = 0, 
         nFeatures.hi = NULL, nFeatures.lo = NULL, 
         dropUnusedDimensions = TRUE, 
         verbose = 2, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Training features (predictive variables). Each column corresponds to a feature and each row to an
observation. 
}
  \item{y}{
The response variable. Can be a single vector or a matrix with arbitrary many columns. Number of rows
(observations) must equal to the number of rows (observations) in x.
}
  \item{xtest}{
Optional test set data. A matrix of the same number of columns (i.e., features) as \code{x}. 
If test set data are not given, only the prediction on training data will be returned. 
}
  \item{classify}{
Should the response be treated as a categorical variable? Classification really only works with two classes.
(The function will run for multiclass problems as well, but the results will be sub-optimal.)
}
  \item{CVfold}{
Optional specification of cross-validation fold. If 0 (the default), no cross-validation is performed.
}
  \item{randomSeed}{
Random seed, used for observation selection for cross-validation. If \code{NULL}, the random generator is
not reset.
}
  \item{assocFnc}{
Function to measure association. Usually a measure of correlation, for example Pearson correlation or
\code{\link{bicor}}. 
}
  \item{assocOptions}{
Character string specifying the options to be passed to the association function.
}
  \item{featureWeightPowers}{
Powers to which to raise the result of \code{assocFnc} to obtain weights. Can be a single number or a vector
of arbitrary length; the returned value will contain one prediction per power.
}
  \item{priorWeights}{
Prior weights for the features. If given, must be either (1) a vector of the same length as the number of
features (columns in \code{x}); (2) a matrix of dimensions length(featureWeightPowers)x(number of features);
or (3) array of dimensions (number of response variables)xlength(featureWeightPowers)x(number of features).
}
  \item{weighByPrediction}{
(Optional) power to downweigh features that are not well predicted between training and test sets. See
details. 
}
  \item{nFeatures.hi}{
Optional restriction of the number of features to use. If given, this many features with the highest association
and lowest association (if \code{nFeatures.lo} is not given) will be used for prediction.
}
  \item{nFeatures.lo}{
Optional restriction of the number of lowest (i.e., most negatively) associated features to use. 
Only used if \code{nFeatures.hi} is also non-NULL.
}
  \item{dropUnusedDimensions}{
Logical: should unused dimensions be dropped from the result?
}
  \item{verbose}{
Integer controling how verbose the diagnostic messages should be. Zero means silent.
}
  \item{indent}{
Indentation for the diagnostic messages. Zero means no indentation, each unit adds two spaces.
}
}
\details{
The predictor calculates the association of each (selected) feature with the response and uses the
association to calculate the weight of the feature as \code{sign(association) *
(association)^featureWeightPower}. Optionally, this weight is multiplied by \code{priorWeights}. Further, a
feature prediction weight can be used to downweigh features that are not well predicted by other features
(see below).

For classification, the (continuous) result of the above calculation is turned into ordinal values
essentially by rounding. 

If features exhibit non-trivial correlations among themselves (such as, for example, in gene expression
data), one can attempt to down-weigh features that do not exhibit the same correlation in the test set.
This is done by using essentially the same predictor to predict _features_ from all other features in the
test data (using the training data to train the feature predictor). Because test features are known, the
prediction accuracy can be evaluated. If a feature is predicted badly (meaning the error in the test set is
much larger than the error in the cross-validation prediction in training data), 
it may mean that its quality in the
training or test data is low (for example, due to excessive noise or outliers). 
Such features can be downweighed using the argument \code{weighByPrediction}. The extra factor is 
min(1, (root mean square prediction error in test set)/(root mean square cross-validation prediction error in
the trainig data)^weighByPrediction), that is it is never bigger than 1.
}
\value{
A list with the following components:
\item{predicted}{The back-substitution prediction on the training data. Normally an array of dimensions
(number of observations) x (number of response variables) x length(featureWeightPowers), but unused
are dropped unless \code{dropUnusedDimensions = FALSE}.}

\item{weightBase}{Absolute value of the associations of each feature with each response.}

\item{variableImportance}{The weight of each feature in the prediction (including the sign).}

\item{predictedTest}{If input \code{xtest} is non-NULL, the predicted test response, in format analogous to \code{predicted} above.}

\item{CVpredicted}{If input \code{CVfold} is non-zero, cross-validation prediction on the training data.}
}
\author{
Peter Langfelder
}
\note{
It makes little practical sense to supply neither \code{xtest} nor \code{CVfold} since the prediction
accuracy on training data will be highly biased.
}

\seealso{
\code{\link{bicor}} for robust correlation that can be used as an association measure
}

\keyword{ misc }
