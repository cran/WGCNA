\name{randomGLMpredictor}
\alias{randomGLMpredictor}
\title{Random generalized linear model predictor }
\description{
An ensemble predictor based on bootstrap aggregation (bagging) of generalized linear models whose covariates are selected using forward stepwise regression according to AIC criteria.
}

\usage{
randomGLMpredictor(
  x, y, xtest = NULL, 
  classify = TRUE,
  nBags = 100,
  replace = TRUE,
  nObsInBag = if (replace) nrow(x) else as.integer(0.632 * nrow(x)),
  nFeaturesInBag = ceiling(ifelse(ncol(x) <= 10, ncol(x), ifelse(ncol(x) <= 300, 
          (1.0276 - 0.00276 * ncol(x)) * ncol(x), ncol(x)/5))),
  nCandidateCovariates=50,
  candidateCorFnc= cor,
  candidateCorOptions = list(method = "pearson", use="p"),
  mandatoryCovariates = NULL,
  randomSeed = 12345,
  verbose =1)
}

\arguments{
  \item{x}{a matrix with rows correspond to observations and columns corresponding to features (covariates).}
  \item{y}{class outcome (factor variable) or quantitative outcome (numeric variable). }
  \item{xtest}{an optional matrix (whose columns correspond to those in x) which contain test (validation) data. The number of rows will typically be different from those in \code{x}.}
  \item{classify}{logical: should the response be treated as 
                  a binary variable (\code{TRUE}) or as a continuous variable (\code{FALSE})?}
  \item{nBags}{number of bags in the ensemble predictor.}
  \item{replace}{logical which deteremines whether the observations for the bag (bootstrap data) are sampled with or without replacement. The function randomly select bagging observations with or without replacement.}
  \item{nObsInBag}{number of observations selected for each bag. Typically, a bootstrap sample (bag) has the same number of observations as in the original observed data (i.e. the rows of \code{x}).}
  \item{nFeaturesInBag}{number of features selected into each bag. Features are randomly selected without replacement.}
  \item{nCandidateCovariates}{top number of features selected with highest absolute correlation with the outcome in individual bag. These features/covariates become the candidates for forward stepwise regression. }
  \item{candidateCorFnc}{the correlation function used to select candidate covariates. Either cor or bicor.}
  \item{candidateCorOptions}{list of arguments to correlation function. If bicor is chosen for class outcome, make sure to include \code{"robustY=F"}.}
  \item{mandatoryCovariates}{indices of features that forced into all regression models across bags. As default, no feature is mandatory.}
  \item{randomSeed}{NULL or integer. The seed for the random number generator. If NULL, the seed will not be set. If non-NULL and the random generator has been initialized prior to the function call, the latter's state is saved and restored upon exit.}
  \item{verbose}{integer which determines the level of verbosity. Zero means silent, higher values make the output progressively more and more verbose.}
}

\details{
The \code{randomGLMpredictor} function requires the R package \code{MASS} since it makes use of the function \code{stepAIC}. Basically, \code{randomGLMpredictor} first selects bootstrapping samples and features randomly for each bag, and then restricts the analysis to features that are highly correlated with the outcome. Prediction in each bag is made based on forward stepwise regression (logistic for binary outcomes, linear for quantitative outcomes). An overall prediction is obtained by averaging results from all bags. Generally, \code{nCandidateCovariates}>100 is not recommended, because the forward selection process is time-consuming. If \code{"nBags=1, replace=F, nObsInBag=nrow(x)"} is used, the function becomes a stepwise generalized linear model predictor without bagging.
}

\value{
  The function returns a list with the following components:
  \item{predictedOOB}{the predicted classification of the input data based on out-of-bag samples. Only for binary outcomes.}
  \item{predictedOOB.cont}{In case of a binary outcome, this is the predicted probability of each outcome specified by \code{y} based on out-of-bag samples. In case of a continous outcome, this is the predicted value based on out-of-bag samples.}
  \item{predictedTest}{if test set is given, the predicted classification for test data. Only for binary outcomes.}
  \item{predictedTest.cont}{if test set is given, the predicted probability of each outcome specified by \code{y} for test data for binary outcomes. In case of a continous outcome, this is the test set predicted value. }
  \item{bagObsIndx}{a matrix with \code{nBags} rows and \code{nObsInBag} columns, giving the indices of observations selected for each bag.}
  \item{datSelectedAsCandidates}{a (0,1) matrix with \code{nBags} rows and columns corresponding to features, indicating which features are selected as candidate regression covariates in each bag.}
  \item{datSelectedByForwardRegression}{a (0,1) matrix with \code{nBags} rows and columns corresponding to features, indicating which features/covariates are selected into the final regression model in each bag.}
  \item{datCoefOfForwardRegression}{a matrix with \code{nBags} rows and columns corresponding to features, giving the final generalized linear model coefficients for features in each bag.}
  \item{timesSelectedByForwardRegression}{a variable importance measure, giving the times each feature is selected into final models in all bags.}
}

\references{
Lin Song, Peter Langfelder, Steve Horvath: Random generalized linear model: a superior ensemble predictor involving few features. BMC Bioinformatics, future.
}

\author{Lin Song
}

\examples{

## binary outcome prediction
# data generation
data(iris)
iris=iris[1:100,]
iris$Species = as.factor(as.character(iris$Species))
set.seed(1)
indx=sample(100, 67, replace=FALSE)
alldat1=iris[indx, ]
alldat2=iris[-indx,]
dat1=alldat1[,-5]
y1=alldat1[,5]
dat2=alldat2[,-5]
y2=alldat2[,5]


# predict with a small number of bags - normally nBags should be at least 100.
RGLM = randomGLMpredictor(dat1, y1, dat2, nCandidateCovariates=ncol(dat1), nBags=30)
y2predict = RGLM$predictedTest
table(y2predict, y2)

}

