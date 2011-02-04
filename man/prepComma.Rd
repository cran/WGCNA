\name{prepComma}
\alias{prepComma}
\title{
Prepend a comma to a non-empty string
}
\description{
Utility function that prepends a comma before the input string if the string is non-empty.
}
\usage{
prepComma(s)
}
\arguments{
  \item{s}{Character string.
}
}
\value{
If \code{s} is non-empty, returns \code{paste(",", s)}, otherwise returns s.
}
\author{
Peter Langfelder
}
\examples{
prepComma("abc");
prepComma("");
}
\keyword{misc}
