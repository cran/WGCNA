\name{colQuantileC}
\alias{colQuantileC}
\alias{rowQuantileC}
\title{ Fast colunm- and row-wise quantile of a matrix. }
\description{
  Fast calculation of column- and row-wise quantiles of a matrix at a single probability. Implemented via compiled
code, it is much faster than the equivalent \code{apply(data, 2, quantile, prob = p)}.
}
\usage{
colQuantileC(data, p)
rowQuantileC(data, p)
}
\arguments{
  \item{data}{ a numerical matrix column-wise quantiles are desired. Missing values are removed.}
  \item{p}{ a single probability at which the quantile is to be calculated. }
}
\details{
At present, only one quantile type is implemented, namely the default type 7 used by R.
}
\value{
  A vector of length equal the number of columns (for \code{colQuantileC}) or rows (for \code{rowQuantileC}) 
  in \code{data} containing the column- or row-wise quantiles.
}

\author{ Peter Langfelder }

\seealso{ 

\code{\link[stats]{quantile}}; 
\code{\link{pquantile}} for another way of calculating quantiles across structured data.
}

\keyword{misc }% __ONLY ONE__ keyword per line
