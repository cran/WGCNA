\name{multiGSub}
\alias{multiGSub}
\alias{multiSub}
\alias{multiGrep}
\alias{multiGrepl}
\title{
Analogs of grep(l) and (g)sub for multiple patterns and relacements
}
\description{
These functions provide convenient pattern finding and substitution for multiple patterns.
}
\usage{
multiGSub(patterns, replacements, x, ...)
multiSub(patterns, replacements, x, ...)
multiGrep(patterns, x, ..., sort = TRUE, invert = FALSE)
multiGrepl(patterns, x, ...)
}
\arguments{
  \item{patterns}{
A character vector of patterns.
}
  \item{replacements}{
A character vector of replacements; must be of the same length as \code{patterns}.
}
  \item{x}{
Character vector of strings in which the pattern finding and replacements should be carried out.
}
  \item{sort}{Logical: should the output indices be sorted in increasing order?}
  \item{invert}{Logical: should the search be inverted and only indices of elements of \code{x} matching none
of the patterns be returned?}
  \item{\dots}{
Other arguments to \code{\link{sub}} or \code{\link{grep}}
}
}
\details{
For each element of \code{x}, patterns are sequentiall searched for and (for \code{multiSub} and
\code{multiGSub} substituted with the corresponding replacement.
}
\value{
\code{multiSub} and \code{multiGSub} return a character vector of the same length as \code{x}, with all
patterns replaces by their replacements in each element of \code{x}. \code{multiSub} replaces each pattern in
each element of \code{x} only once, \code{multiGSub} as many times as the pattern is found.

\code{multiGrep} returns the indices of those elements in \code{x} in which at least one of \code{patterns}
was found, or, if \code{invert} is TRUE, the indices of elements in which none of the patterns were found. 

\code{multiGrepl} returns a logical vector of the same length as \code{x}, with \code{TRUE} is any of the
patterns matched the element of \code{x}, and \code{FALSE} otherwise. 
}
\author{
Peter Langfelder
}
\seealso{
The workhorse functions \code{\link{sub}}, \code{\link{gsub}}, \code{\link{grep}} and \code{\link{grepl}}.
}
\keyword{misc}% __ONLY ONE__ keyword per line
