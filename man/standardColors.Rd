\name{standardColors}
\alias{standardColors}
\title{Colors this library uses for labeling modules.} 
\description{
  Returns the vector of color names in the order they are assigned by other functions in this library.
}
\usage{
standardColors(n = NULL)
}
\arguments{
  \item{n}{Number of colors requested. If \code{NULL}, all (approx. 450) colors will be returned. Any
other invalid argument such as less than one or more than maximum (\code{length(standardColors())}) will
trigger an error. }
}
\value{
  A vector of character color names of the requested length.
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\examples{
standardColors(10);
}
\keyword{color}
\keyword{misc}
