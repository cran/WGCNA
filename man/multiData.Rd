\name{multiData}
\alias{multiData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a multiData structure.
}
\description{
This function creates a multiData structure by storing its input arguments as the 'data' components.
}
\usage{
multiData(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
  Arguments to be stored in the multiData structure.
}
}
\details{
A multiData structure is intended to store (the same type of) data for multiple, possibly independent,
realizations
(for example, expression data for several independent experiments). It is a list where
each component corresponds to an (independent) data set. Each component is in turn a list that can hold
various types of information but must have a \code{data} component. In a "strict" multiData structure, the
\code{data} components are required to each be a matrix or a data frame and have the same number of
columns. In a "loose" multiData structure, the \code{data} components can be anything (but for most
purposes should be of comparable type and content).

}
\value{
  The resulting multiData structure. 
}
\author{
Peter Langfelder
}

\seealso{
  \code{\link{multiData2list}} for converting a multiData structure to a list;
  \code{\link{list2multiData}} for an alternative way of creating a multiData structure;
  \code{\link{mtd.apply}, \link{mtd.applyToSubset}, \link{mtd.mapply}} for ways of applying a function to
each component of a multiData structure.
}
\examples{
data1 = matrix(rnorm(100), 20, 5);
data2 = matrix(rnorm(50), 10, 5);

md = multiData(Set1 = data1, Set2 = data2);

checkSets(md)
}
\keyword{misc }% __ONLY ONE__ keyword per line
