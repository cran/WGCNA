\name{nSets}
\alias{nSets}
\title{ Number of sets in a multi-set variable
}
\description{
A convenience function that returns the number of sets in a multi-set variable.
}
\usage{
nSets(multiData, ...)
}
\arguments{
  \item{multiData}{
vector of lists; in each list there must be a component named \code{data} whose content
is a matrix or dataframe or array of dimension 2. }
  \item{\dots}{
Other arguments to function \code{\link{checkSets}}. 
}
}
\value{
  A single integer that equals the number of sets given in the input \code{multiData}.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{checkSets}}
}
\keyword{misc}
