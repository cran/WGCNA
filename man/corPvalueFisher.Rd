\name{corPvalueFisher}
\alias{corPvalueFisher}
\title{ Fisher's asymptotic p-value for correlation}
\description{
  Calculates Fisher's asymptotic p-value for given correlations. 
}
\usage{
corPvalueFisher(cor, nSamples, twoSided = TRUE)
}
\arguments{
  \item{cor}{ A vector of correlation values whose corresponding p-values are to be calculated }
  \item{nSamples}{ Number of samples from which the correlations were calculated }
  \item{twoSided}{ logical: should the calculated p-values be two sided? }
}
\value{
  A vector of p-values of the same length as the input correlations. 
}
\author{ Steve Horvath and Peter Langfelder }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
