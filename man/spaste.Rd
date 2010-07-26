\name{spaste}
\alias{spaste}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Space-less paste
}
\description{
  A convenient wrapper for the \code{\link{paste}} function with \code{sep=""}.
}
\usage{
spaste(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
standard arguments to function \code{\link{paste}} except \code{sep}.
}
}
\value{
  The result of the corresponding \code{\link{paste}}.
}
\author{
Peter Langfelder
}
\note{
Do not use the \code{sep} argument. Using will lead to an error.
}

\seealso{
\code{\link{paste}}
}
\examples{
  a = 1;
  paste("a=", a);
  spaste("a=", a);
}
\keyword{misc}
