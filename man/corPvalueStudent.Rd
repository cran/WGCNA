\name{corPvalueStudent}
\alias{corPvalueStudent}
\title{Student asymptotic p-value for correlation}
\description{
  Calculates Student asymptotic p-value for given correlations. 
}
\usage{
corPvalueStudent(cor, nSamples)
}
\arguments{
  \item{cor}{ A vector of correlation values whose corresponding p-values are to be calculated }
  \item{nSamples}{ Number of samples from which the correlations were calculated }
}

\value{
 A vector of p-values of the same length as the input correlations. 
}
\author{ Steve Horvath and Peter Langfelder }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
