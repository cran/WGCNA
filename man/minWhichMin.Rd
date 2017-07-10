\name{minWhichMin}
\alias{minWhichMin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast joint calculation of row- or column-wise minima and indices of minimum elements
}
\description{
Fast joint calculation of row- or column-wise minima and indices of minimum elements. Missing data are
removed.
}
\usage{
minWhichMin(x, byRow = FALSE, dims = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric matrix or array.
}
  \item{byRow}{
Logical: should the minima and indices be found for columns (\code{FALSE}) or rows (\code{TRUE})?
}
  \item{dims}{
Specifies dimensions for which to find the minima and indices. For \code{byRow = FALSE}, they are calculated for
dimensions \code{dims+1} to \code{n=length(dim(x))}; for For \code{byRow = TRUE}, they are calculated for dimensions
1,...,\code{dims}.
}
}
\value{
A list with two components, \code{min} and \code{which}; each is a vector or array with dimensions 

\code{dim(x)[(dims+1):n]} (with \code{n=length(dim(x))}) if \code{byRow = FALSE}, and

\code{dim(x)[1:dims]} if \code{byRow = TRUE}.
}
\author{
Peter Langfelder
}
\keyword{stats}% __ONLY ONE__ keyword per line
