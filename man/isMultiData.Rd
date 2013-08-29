\name{isMultiData}
\alias{isMultiData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine whether the supplied object is a valid multiData structure
}
\description{
Attempts to determine whether the supplied object is a valid multiData structure (see Details).
}
\usage{
isMultiData(x, strict = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object.
}
  \item{strict}{Logical: should the structure of multiData be checked for "strict" compliance?}
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

This function checks whether the supplied \code{x} is a multiData structure in the "strict" (when
\code{strict = TRUE} or "loose" \code{strict = FALSE} sense.

}
\value{
Logical: \code{TRUE} if the input \code{x} is a multiData structure, \code{FALSE} otherwise.
}
\author{
Peter Langfelder
}

\seealso{
Other multiData handling functions whose names start with \code{mtd.}
}
\keyword{ misc}% __ONLY ONE__ keyword per line
