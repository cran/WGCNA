\name{relativeCorPredictionSuccess}
\alias{relativeCorPredictionSuccess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare prediction success }
\description{
  Compare prediction success of several gene screening methods.
}
\usage{
relativeCorPredictionSuccess(
  corPredictionNew, 
  corPredictionStandard, 
  corTestSet, 
  topNumber = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corPredictionNew}{ Matrix of predictor statistics }
  \item{corPredictionStandard}{ Reference presdictor statistics}
  \item{corTestSet}{ Correlations of predictor variables with trait in test set}
  \item{topNumber}{ A vector giving the numbers of top genes to consider }
}
\value{
  Data frame with components
  \item{topNumber}{copy of the input \code{topNumber}}
  \item{kruskalp}{Kruskal-Wallis p-values}
}
\author{ Steve Horvath }

\seealso{ \code{\link{corPredictionSuccess}} }
\keyword{ misc }
