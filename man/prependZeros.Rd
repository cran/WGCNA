\name{prependZeros}
\alias{prependZeros}
\alias{prependZeros.int}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Pad numbers with leading zeros to specified total width
}
\description{
  These functions pad the specified numbers with zeros to a specified total width.
}
\usage{
prependZeros(x, width = max(nchar(x)))
prependZeros.int(x, width = max(nchar(as.integer(x))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     Vector of numbers to be padded. For \code{prependZeros}, the vector may be real (non-integer) or even character (and
not necessarily representing numbers). For \code{prependZeros}, the vector must be numeric and non-integers get rounded
down to the nearest integer.
}
  \item{width}{
     Width to pad the numbers to.
}
}
\details{
  The \code{prependZeros.int} version works better with numbers such as 100000 which may get converted to character as 1e5 and
hence be incorrectly padded in the \code{prependZeros} function. On the flip side, prependZeros works also for non-integer
inputs.
}
\value{
  Character vector with the 0-padded numbers.
}
\author{
  Peter Langfelder
}
\examples{
prependZeros(1:10)
prependZeros(1:10, 4)
# more exotic examples
prependZeros(c(1, 100000), width = 6) ### Produces incorrect output
prependZeros.int(c(1, 100000))  ### Correct output
prependZeros(c("a", "b", "aa")) ### pads the shorter strings using zeros.
}
\keyword{ misc }
