\name{colQuantileC}
\alias{colQuantileC}
\title{ Fast colunm-wise quantile of a matrix. }
\description{
  Fast calculation of column-wise quantiles of a matrix at a single probability. Implemented via compiled
code, it is much faster than the equivalent \code{apply(data, 2, quantile, prob = p)}.
}
\usage{
colQuantileC(data, p)
}
\arguments{
  \item{data}{ a numerical matrix column-wise quantiles are desired. Missing values are corrently not
allowed. }
  \item{p}{ a single probability at which the quantile is to be calculated. }
}
\value{
  A vector of length equal the number of columns in \code{data} containing the column-wise quantiles.
}

\author{ Peter Langfelder }

\seealso{ \code{\link[stats]{quantile}} }

\keyword{misc }% __ONLY ONE__ keyword per line
