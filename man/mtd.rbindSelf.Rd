\name{mtd.rbindSelf}
\alias{mtd.rbindSelf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turn a multiData structure into a single matrix or data frame.
}
\description{
This function "rbinds" the \code{data} components of all sets in the input into a single matrix or data
frame.
}
\usage{
mtd.rbindSelf(multiData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiData}{
A multiData structure.
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

This function requires a "strict" multiData structure.

}
\value{
A single matrix or data frame containing the "rbinded" result.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{rbind}} for various subtleties of the row binding operation.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
