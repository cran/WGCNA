\name{nearestCentroidPredictor}
\alias{nearestCentroidPredictor}
\title{
Nearest centroid predictor
}
\description{
Nearest centroid predictor for binary (i.e., two-outcome) data. Implements a whole host of options and
improvements such as accounting for within-class heterogeneity using sample networks, 
various ways of feature selection and weighing etc.
}
\usage{
nearestCentroidPredictor(

    # Input training and test data
    x, y, 
    xtest = NULL, 

    # Feature weights and selection criteria
    featureSignificance = NULL, 
    assocFnc = "cor", assocOptions = "use = 'p'", 
    assocCut.hi = NULL, assocCut.lo = NULL, 
    nFeatures.hi = 10, nFeatures.lo = 10,
    weighFeaturesByAssociation = 0, 
    scaleFeatureMean = TRUE, scaleFeatureVar = TRUE, 

    # Predictor options 
    centroidMethod = c("mean", "eigensample"), 
    simFnc = "cor", simOptions = "use = 'p'", 
    useQuantile = NULL, 
    sampleWeights = NULL, 
    weighSimByPrediction = 0, 

    # What should be returned
    CVfold = 0, returnFactor = FALSE, 

    # General options
    randomSeed = 12345, 
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
  \item{featureSignificance}{
Optional vector of feature significance for the response variable. If given, it is used for feature
selection (see details). Should preferably be signed, that is features can have high negative significance.
}
  \item{assocFnc}{
Character string specifying the association function. The association function should behave roughly as
\code{link{cor}} in that it takes two arguments
(a matrix and a vector) plus options 
and returns the vector of associations between the columns of the matrix and the
vector. The associations may be signed (i.e., negative or positive). 
}
  \item{assocOptions}{
Character string specifying options to the association function.
}
  \item{assocCut.hi}{
Association (or featureSignificance) threshold for including features in the predictor. Features with
associtation higher than \code{assocCut.hi} will be included. If not given, the threshold method will not be
used; instead, a fixed number of features will be included as specified by \code{nFeatures.hi}
and \code{nFeatures.lo}.
}
  \item{assocCut.lo}{
Association (or featureSignificance) threshold for including features in the predictor. Features with 
associtation lower than \code{assocCut.lo} will be included. If not given, defaults to \code{-assocCut.hi}.
If \code{assocCut.hi} is \code{NULL}, the threshold method will not be
used; instead, a fixed number of features will be included as specified by \code{nFeatures.hi} and
\code{nFeatures.lo}. 
}
  \item{nFeatures.hi}{
Number of highest-associated features (or features with highest \code{featureSignificance}) to include in the
predictor. Only used if \code{assocCut.hi} is \code{NULL}.
}
  \item{nFeatures.lo}{
Number of lowest-associated features (or features with highest \code{featureSignificance}) to include in
the predictor. Only used if \code{assocCut.hi} is \code{NULL}.
}
  \item{weighFeaturesByAssociation}{
(Optional) power to downweigh features that are less associated with the response. See details.
}
  \item{scaleFeatureMean}{
Logical: should the training features be scaled to mean zero? Unless there are good reasons not to scale,
the features should be scaled.
}
  \item{scaleFeatureVar}{
Logical: should the training features be scaled to unit variance? Again, unless there are good reasons not
to scale, the features should be scaled.
}
  \item{centroidMethod}{
One of \code{"mean"} and \code{"eigensample"}, specifies how the centroid should be calculated.
\code{"mean"} takes the mean across all samples (or all samples within a sample module, if sample networks
are used), whereas \code{"eigensample"} calculates the first principal component of the feature matrix and
uses that as the centroid.
}
  \item{simFnc}{
Character string giving the similarity function for measuring the similarity between test samples and
centroids. This function should
behave roughly like the function \code{\link{cor}} in that it takes two arguments (\code{x}, \code{y})
and calculates the pair-wise similarities between columns of \code{x} and \code{y}. For convenience, the
value \code{"dist"} is treated specially: the Euclidean distance between the columns of \code{x} and
\code{y} is calculated and its negative is returned (so that smallest distance corresponds to highest
similarity). Since values of this function are only used for ranking centroids, its values are not
restricted to be positive or within certain bounds.
}
  \item{simOptions}{
Character string specifying the options to the similarity function.
}
  \item{useQuantile}{
If non-NULL, the "nearest quantiloid" will be used instead of the nearest centroid. See details.
}
  \item{sampleWeights}{
Optional specification of sample weights. Useful for example if one wants to explore boosting.
}
  \item{weighSimByPrediction}{
(Optional) power to downweigh features that are not well predicted between training and test sets. See
details.
}
  \item{CVfold}{
Non-negative integer specifying cross-validation. Zero means no cross-validation will be performed. values
above zero specify the number of samples to be considered test data for each step of cross-validation.
}
  \item{returnFactor}{
Logical: should a factor be returned?
}
  \item{randomSeed}{
Integere specifying the seed for the random number generator. If \code{NULL}, the seed will not be set. See
\code{\link{set.seed}}.
}
  \item{verbose}{
Integer controling how verbose the diagnostic messages should be. Zero means silent.
}
  \item{indent}{
Indentation for the diagnostic messages. Zero means no indentation, each unit adds two spaces.
}
}
\details{

Nearest centroid predictor works by forming a representative profile (centroid) 
across features for each class from
the training data, then assigning each test sample to the class of the nearest representative profile. The
representative profile can be formed either as mean or as athe first principal component ("eigensample";
this choice is governed by the option \code{centroidMethod}).

When the number of features is large and only a small fraction is likely to be associated with the outcome, 
feature selection can be used to restrict the features that actually enter the centroid. Feature selection
can be based either on their association with the outcome
calculated from the training data using \code{assocFnc}, or on user-supplied feature significance (e.g.,
derived from literature, argument
\code{featureSignificance}). In either case, features can be selected by high and low association tresholds 
or by taking a fixed number of highest- and lowest-associated features. 

As an alternative to centroids, the predictor can also assign test samples based on a given quantile of the
distances from the training samples in each class (argument \code{useQuantile}). This may be advantageous if
the samples in each class form irregular clusters. Note that setting \code{useQuantile=0} (i.e., using
minimum distance in each class) essentially gives a nearest neighbor predictor: each test sample will be
assigned to the class of its nearest training neighbor.

If features exhibit non-trivial correlations among themselves (such as, for example, in gene expression
data), one can attempt to down-weigh features that do not exhibit the same correlation in the test set.
This is done by using essentially the same predictor to predict _features_ from all other features in the
test data (using the training data to train the feature predictor). Because test features are known, the
prediction accuracy can be evaluated. If a feature is predicted badly (meaning the error in the test set is
much larger than the error in the cross-validation prediction in training data),
it may mean that its quality in the
training or test data is low (for example, due to excessive noise or outliers).
Such features can be downweighed using the argument \code{weighByPrediction}. The extra factor is
min(1, (root mean square prediction error in test set)/(root mean square cross-validation prediction error
in
the trainig data)^weighByPrediction), that is it is never bigger than 1.

Unless the features' mean and variance can be ascribed clear meaning, the (training) features should be
scaled to mean 0 and variance 1 before the centroids are formed. 

The function implements a basic option for removal of spurious effects in the training and test data, by
removng a fixed number of leading principal components from the features. This sometimes leads to better
prediction accuracy but should be used with caution.

If samples within each class are heterogenous, a single centroid may not represent each class well. This
function can deal with within-class heterogeneity by clustering samples (separately in each class), then
using a one representative (mean, eigensample) or quantile for each cluster in each class to assign test
samples. Various similarity measures, specified by \code{adjFnc}, can be used to construct the sample network
adjacency. Similarly, the user can specify a clustering function using \code{clusteringFnc}. The
requirements on the clustering function are described in a separate section below.

}
\value{
A list with the following components:

\item{predicted}{The back-substitution prediction in the training set.}

\item{predictedTest}{Prediction in the test set.}

\item{featureSignificance}{A vector of feature significance calculated by \code{assocFnc} or a copy of the
input \code{featureSignificance} if the latter is non-NULL.}

\item{selectedFeatures}{A vector giving the indices of the features that were selected for the predictor.}

\item{centroidProfile}{The representative profiles of each class (or cluster). Only returned in
\code{useQuntile} is \code{NULL}. }

\item{testSample2centroidSimilarities}{A matrix of calculated similarities between the test samples and
class/cluster centroids.}

\item{featureValidationWeights}{A vector of validation weights (see Details) for the selected features. If
\code{weighFeaturesByValidation} is 0, a unit vector is used and returned.}

\item{CVpredicted}{Cross-validation prediction on the training data. Present only if \code{CVfold} is
non-zero.}

\item{sampleClusterLabels}{A list with two components (one per class). Each component is a vector of sample
cluster labels for samples in the class.}

}
\author{
Peter Langfelder
}
\seealso{
\code{\link{votingLinearPredictor}}
}
\keyword{misc} 
