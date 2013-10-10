\name{mtd.mapply}
\alias{mtd.mapply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply a function to elements of given multiData structures.
}
\description{
Inspired by \code{\link{mapply}}, this function applies a given function to each \code{data} component in
the input multiData arguments, and optionally simplify the result to an array if possible.

}
\usage{
mtd.mapply(FUN, ..., MoreArgs = NULL, 
           mdmaSimplify = FALSE, 
           mdma.doCollectGarbage = FALSE, 
           mdma.argIsMultiData = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FUN}{Function to be applied. }

  \item{\dots}{
   Arguments to be vectorized over. These can be multiData structures or simple vectors (e.g., lists). 
}
  \item{MoreArgs}{
   A named list that specifies the scalar arguments (if any) to \code{FUN}.
}
  \item{mdmaSimplify}{
   Logical: should simplification of the result to an array be attempted? The simplification is fragile and
can produce unexpected errors; use the default \code{FALSE} if that happens.
}
  \item{mdma.doCollectGarbage}{
   Should garbage collection be forced after each application of \code{FUN}?
}
  \item{mdma.argIsMultiData}{
   Optional specification whether arguments are multiData structures. A logical vector where each component
corresponds to one entry of \code{...}. If not given, multiData status will be determined using
\code{\link{isMultiData}} with argument \code{strict=FALSE}.
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

This function applies the function \code{FUN} to each \code{data} component of those arguments in
\code{...} that are multiData structures in the "loose" sense, 
and to each component of those arguments in \code{...} that are
not multiData structures.
}
\value{
A multiData structure containing (as the \code{data} components) the results of \code{FUN}. If
simplification is successful, an array instead.
}
\author{
Peter Langfelder
}
\seealso{
\code{multiData.apply} for application of a function to a single multiData structure.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
