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
mtd.mapply(

  # What to do
  FUN, ..., MoreArgs = NULL, 

  # How to interpret the input
  mdma.argIsMultiData = NULL,

  # Copy previously known results?
  mdmaExistingResults = NULL, mdmaUpdateIndex = NULL,

  # How to format output
  mdmaSimplify = FALSE, 
  returnList = FALSE,
  
  # Options controlling internal behaviour
  mdma.doCollectGarbage = FALSE, 
  mdmaVerbose = 0, mdmaIndent = 0)

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

  \item{mdma.argIsMultiData}{
   Optional specification whether arguments are multiData structures. A logical vector where each component
corresponds to one entry of \code{...}. If not given, multiData status will be determined using
\code{\link{isMultiData}} with argument \code{strict=FALSE}.
}

  \item{mdmaExistingResults}{Optional list that contains previously calculated results. This can be useful
if only a few sets in \code{multiData} have changed and recalculating the unchanged ones is computationally
expensive. If not given, all calculations will be performed. If given, components of this list are copied
into the output. See \code{mdmUpdateIndex} for which components are re-calculated by default. }

  \item{mdmaUpdateIndex}{Optional specification of which sets in \code{multiData} the calculation should 
actually be carried out. This argument has an effect only if \code{mdmaExistingResults} is non-NULL. If the
length of \code{mdmaExistingResults} (call the length `k') 
is less than the number of sets in \code{multiData}, the function
assumes that the existing results correspond to the first `k' sets in \code{multiData} and the rest of the
sets are automatically calculated, irrespective of the setting of \code{mdmaUpdateIndex}. The argument
\code{mdmaUpdateIndex} can be used to specify re-calculation of some (or all) of the results that already
exist in \code{mdmaExistingResults}. }

  \item{mdmaSimplify}{
   Logical: should simplification of the result to an array be attempted? The simplification is fragile and
can produce unexpected errors; use the default \code{FALSE} if that happens.
}

  \item{returnList}{Logical: should the result be turned into a list (rather than a multiData structure)?
Note that this is incompatible with simplification: if \code{mdaSimplify} is \code{TRUE}, this argument is
ignored.}

  \item{mdma.doCollectGarbage}{
   Should garbage collection be forced after each application of \code{FUN}?
}
  \item{mdmaVerbose}{Integer specifying whether progress diagnistics should be printed out. Zero means
silent, increasing values will lead to more diagnostic messages.}

  \item{mdmaIndent}{Integer specifying the indentation of the printed progress messages. Each unit equals
two spaces.}

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
\code{\link{multiData}} to create a multiData structure;

\code{multiData.apply} for application of a function to a single multiData structure.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
