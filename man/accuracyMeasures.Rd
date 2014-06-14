\name{accuracyMeasures}
\alias{accuracyMeasures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accuracy measures for a 2x2 confusion matrix or for vectors of predicted and observed values.
}
\description{
The function calculates various prediction accuracy statistics for predictions of binary or quantitative
(continuous) responses. For binary classification, the function calculates 
the error rate, accuracy, sensitivity, specificity, positive predictive value, and
other accuracy measures. For quantitative prediction, the function calculates correlation, R-squared, error
measures, and the C-index.
}
\usage{
accuracyMeasures(
  predicted, 
  observed = NULL, 
  type = c("auto", "binary", "quantitative"),
  levels = if (isTRUE(all.equal(dim(predicted), c(2,2)))) colnames(predicted)
            else if (is.factor(predicted))
              sort(unique(c(as.character(predicted), as.character(observed))))
            else sort(unique(c(observed, predicted))),
  negativeLevel = levels[2], 
  positiveLevel = levels[1] )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predicted}{ either a a 2x2 confusion matrix (table) whose entries contain non-negative
integers, or a vector of predicted values. Predicted values can be binary or quantitative (see \code{type}
below). If a 2x2 matrix is given, it must have valid column and row names that specify the levels of the
predicted and observed variables whose counts the matrix is giving (e.g., the
function \code{\link{table}} sets the names appropriately.) If it is a 2x2 table and the table 
contains non-negative real (non-integer) numbers the function outputs a warning.
}
  \item{observed}{ if \code{predicted} is a vector of predicted values, this (\code{observed}) must be a
vector of the same length giving the "gold standard" (or observed) values. Ignored if \code{predicted} is a 2x2
table. }

  \item{type}{ character string specifying the type of the prediction problem (i.e., values in the 
\code{predicted} and \code{observed} vectors). The
default \code{"auto"} decides type automatically: 
if \code{predicted} is a 2x2 table or if the number of unique values in
the concatenation of \code{predicted} and \code{observed} is 2, the prediction problem (type) is assumed to
be binary, otherwise it is assumed to be quantitative. Inconsistent specification (for example, when
\code{predicted} is a 2x2 matrix and \code{type} is \code{"quantitative"}) trigger errors. } 

  \item{levels}{ a 2-element vector specifying the two levels of binary variables. Only used if \code{type} is
\code{"binary"} (or \code{"auto"} that results in the binary type). Defaults to either the column names of
the confusion matrix (if the matrix is specified) or to the sorted unique values of \code{observed} and
\code{opredicted}. }

  \item{negativeLevel}{ the binary value (level) that corresponds to the negative outcome. Note that the
default is the second of the sorted levels (for example, if levels are 1,2, the default negative level is
2). Only used if \code{type} is
\code{"binary"} (or \code{"auto"} that results in the binary type).}

   \item{positiveLevel}{ the binary value (level) that corresponds to the positive outcome. Note that the 
default is the second of the sorted levels (for example, if levels are 1,2, the default negative level is
2). Only used if \code{type} is
\code{"binary"} (or \code{"auto"} that results in the binary type).}

}
\details{ 


The rows of the 2x2 table tab must correspond to a test (or predicted) outcome and the columns to a true
outcome ("gold standard"). A table that relates a predicted outcome to a true test outcome is also known as
confusion matrix. Warning: To correctly calculate sensitivity and specificity, the positive and negative
outcome must be properly specified so they can be matched to the appropriate rows and columns in the
confusion table. 

Interchanging the negative and positive levels swaps the estimates of the sensitivity and specificity 
but has no effect on the error rate or
accuracy. Specifically, denote by \code{pos} the index of the positive level in the confusion table, and by
\code{neg} th eindex of the negative level in the confusion table. 
The function then defines number of true positives=TP=tab[pos, pos], no.false positives
=FP=tab[pos, neg], no.false negatives=FN=tab[neg, pos], no.true negatives=TN=tab[neg, neg]. 
Then Specificity= TN/(FP+TN)
Sensitivity= TP/(TP+FN) NegativePredictiveValue= TN/(FN + TN) PositivePredictiveValue= TP/(TP + FP)
FalsePositiveRate = 1-Specificity FalseNegativeRate = 1-Sensitivity Power = Sensitivity
LikelihoodRatioPositive = Sensitivity / (1-Specificity) LikelihoodRatioNegative =
(1-Sensitivity)/Specificity. The naive error rate is the error rate of a constant (naive) predictor that
assigns the same outcome to all samples. The prediction of the naive predictor equals the most frequenly
observed outcome. Example: Assume you want to predict disease status and 70 percent of the observed samples
have the disease. Then the naive predictor has an error rate of 30 percent (since it only misclassifies 30
percent of the healthy individuals). }
\value{
Data frame with two columns: 
\item{Measure}{this column contais character strings that specify name of the accuracy measure.}
\item{Value}{this column contains the numeric estimates of the corresponding accuracy measures.}
}
\references{
http://en.wikipedia.org/wiki/Sensitivity_and_specificity 
}
\author{
Steve Horvath and Peter Langfelder
}
\examples{
m=100
trueOutcome=sample( c(1,2),m,replace=TRUE)
predictedOutcome=trueOutcome
# now we noise half of the entries of the predicted outcome
predictedOutcome[ 1:(m/2)] =sample(predictedOutcome[ 1:(m/2)] )
tab=table(predictedOutcome, trueOutcome) 
accuracyMeasures(tab)

# Same result:
accuracyMeasures(predictedOutcome, trueOutcome)

}
\keyword{ misc }
