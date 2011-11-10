\name{accuracyMeasures}
\alias{accuracyMeasures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accuracy measures for a 2x2 confusion matrix.
}
\description{
The function calculates the error rate, accuracy, sensitivity, specificity, positive predictive value, and other accuracy measures for a 2x2 confusion matrix. The input \code{tab} must be a 2x2 dimensional matrix which contains count data. 
}
\usage{
accuracyMeasures(tab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{ should be a 2x2 dimensional matrix (i.e. table) whose entries contain non-negative integers. The function also works if the table contains non-negative real numbers but outputs a warning if non-integers are detected. 
}
}
\details{
The rows of the 2x2 table tab may correspond to a test (or predicted) outcome and the columns to a true outcome ("gold standard"). A table that relates a predicted outcome to a true test outcome is also known as confusion matrix. Warning: When it comes to estimating the sensitivity and specificity, the order of rows and columns matters. The first row corresponds to a *positive* outcome while the second row corresponds to a negative outcome. Same holds for the columns. Reversing the order of the rows (and accordingly the order of the columns) swaps the estimates of the sensitivity and specificity but has no effect on the error rate or accuracy. Specifically, the funcdtion assumes: number of true positives=TP=tab[1,1], no.false positives =FP=tab[1,2], no.false negatives=FN=tab[2,1], no.true negatives=TN=tab[2,2]. Then Specificity= TN/(FP+TN) Sensitivity= TP/(TP+FN) NegativePredictiveValue= TN/(FN + TN) PositivePredictiveValue= TP/(TP + FP) FalsePositiveRate = 1-Specificity FalseNegativeRate = 1-Sensitivity Power = Sensitivity LikelihoodRatioPositive = Sensitivity / (1-Specificity) LikelihoodRatioNegative = (1-Sensitivity)/Specificity. The naive error rate is the error rate of a constant (naive) predictor that assigns the same outcome to all samples. The prediction of the naive predictor equals the most frequenly observed outcome. Example: Assume you want to predict disease status and 70 percent of the observed samples have the disease. Then the naive predictor has an error rate of 30 percent (since it only misclassifies 30 percent of the healthy individuals). }
\value{
Data frame with two columns: 
\item{Measure}{this column contais character strings that specify name of the accuracy measure.}
\item{Value}{this column contains the numeric estimates of the corresponding accuracy measures.}
}
\references{
http://en.wikipedia.org/wiki/Sensitivity_and_specificity 
}
\author{
Steve Horvath
}
\examples{
m=100
trueOutcome=sample( c(1,2),m,replace=TRUE)
predictedOutcome=trueOutcome
# now we noise half of the entries of the predicted outcome
predictedOutcome[ 1:(m/2)] =sample(predictedOutcome[ 1:(m/2)] )
tab=table(predictedOutcome, trueOutcome) 
accuracyMeasures(tab)

}
\keyword{ misc }
