\name{mtd.apply}
\alias{mtd.apply}
\alias{mtd.applyToSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply a function to each set in a multiData structure.
}
\description{
Inspired by \code{\link{lapply}}, these functions apply a given function to each \code{data} component in
the input \code{multiData} structure, and optionally simplify the result to an array if possible.
}
\usage{
mtd.apply(multiData, FUN, ..., 
          mdaSimplify = FALSE, mdaCopyNonData = FALSE)

mtd.applyToSubset(multiData, FUN, ..., 
                  mdaRowIndex = NULL, mdaColIndex = NULL, 
                  mdaSimplify = FALSE, mdaCopyNonData = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiData}{
A multiData structure to apply the function over
}
  \item{FUN}{Function to be applied. }

  \item{\dots}{
Other arguments to the function \code{FUN}.
}
  \item{mdaRowIndex}{If given, must be a list of the same length as \code{multiData}. Each element must be
a logical or numeric vector that specifies rows in each \code{data} component 
               to select before applying the function.}
  \item{mdaColIndex}{A logical or numeric vector that specifies columns in each \code{data} component 
               to select before applying the function.}
  \item{mdaSimplify}{
Logical: should the result be simplified to an array, if possible? Note that this may lead to errors; if
so, disable simplification.
}
  \item{mdaCopyNonData}{Logical: should non-data components of \code{multiData} be copied into the output?
Note that the copying is incompatible with simplification; enabling both will trigger an error.}
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

\code{mtd.apply} works on any "loose" multiData structure; \code{mtd.applyToSubset} assumes (and checks
for) a "strict" multiData structure.

}
\value{
A multiData structure containing the results of the supplied function on each \code{data} component in the
input multiData structure. Other components are simply copied.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{mtd.applyToSubset}} for applying a function to a subset of a multiData structure;
\code{\link{mtd.mapply}} for vectorizing over several arguments.

}
\keyword{ misc }% __ONLY ONE__ keyword per line
