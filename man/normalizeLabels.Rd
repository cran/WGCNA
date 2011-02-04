\name{normalizeLabels}
\alias{normalizeLabels}
\title{Transform numerical labels into normal order. }
\description{
  Transforms numerical labels into normal order, that is the largest group will be labeled 1, next
largest 2 etc. Label 0 is optionally preserved.
}
\usage{
normalizeLabels(labels, keepZero = TRUE)
}
\arguments{
  \item{labels}{Numerical labels.}
  \item{keepZero}{If \code{TRUE} (the default), labels 0 are preserved.}
}
\value{
  A vector of the same length as input, containing the normalized labels.
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
