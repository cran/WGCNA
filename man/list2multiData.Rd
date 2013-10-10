\name{list2multiData}
\alias{list2multiData}
\alias{multiData2list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a list to a multiData structure and vice-versa.
}
\description{
\code{list2multiData} converts a list to a multiData structure; \code{multiData2list} does the inverse.
}
\usage{
list2multiData(data)
multiData2list(multiData)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A list to be converted to a multiData structure.
}
  \item{multiData}{ A multiData structure to be converted to a list. }
}
\details{
A multiData structure is a vector of lists (one list for each set) where each list has a
component \code{data} containing some useful information. 
}
\value{
For \code{list2multiData}, a multiData structure; for \code{multiData2list}, the corresponding list.
}
\author{
Peter Langfelder
}

\keyword{misc}% __ONLY ONE__ keyword per line
