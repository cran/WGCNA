\name{prependZeros}
\alias{prependZeros}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Pad numbers with leading zeros to specified total width
}
\description{
  This function pads the specified numbers with zeros to a specified total width.
}
\usage{
prependZeros(x, width = max(nchar(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     Vector of numbers to be padded.
}
  \item{width}{
     Width to pad the numbers to.
}
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
}
\keyword{ misc }
