\name{corPredictionSuccess}
\alias{corPredictionSuccess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Qunatification of success of gene screening }
\description{
  This function calculates the success of gene screening.
}
\usage{
corPredictionSuccess(corPrediction, corTestSet, topNumber = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corPrediction}{ a vector or a matrix of prediction statistics }
  \item{corTestSet}{ correlation or other statistics on test set }
  \item{topNumber}{ a vector of the number of top genes to consider}
}
\details{
  For each column in \code{corPrediction}, the function evaluates the mean \code{corTestSet} for the number
of top genes (ranked by the column in \code{corPrediction}) given in \code{topNumber}. The higher the mean
\code{corTestSet} (for positive \code{corPrediction}) or negative (for negative \code{corPrediction}), the
more successful the prediction.
}
\value{
  \item{meancorTestSetOverall }{ difference of \code{meancorTestSetPositive} and
\code{meancorTestSetNegative} below }
  \item{meancorTestSetPositive}{ mean \code{corTestSet} on top genes with positive \code{corPrediction} }
  \item{meancorTestSetNegative}{ mean \code{corTestSet} on top genes with negative \code{corPrediction} }
  ...
}
\author{ Steve Horvath }

\seealso{ \code{\link{relativeCorPredictionSuccess}}}
\keyword{ misc}
