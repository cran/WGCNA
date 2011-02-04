\name{matchLabels}
\alias{matchLabels}
\title{ Relabel module labels to best match the given reference labels }
\description{
  Given a \code{source} and  \code{reference} vectors of module labels, the function produces a module
labeling that is equivalent to \code{source}, but individual modules are re-labeled so that modules with
significant overlap in \code{source} and \code{reference} have the same labels.
}
\usage{
matchLabels(source, reference, pThreshold = 5e-2)
}
\arguments{
  \item{source}{ a vector or a matrix of reference labels. The labels may be numeric or character. }
  \item{reference}{ a vector of reference labels. }
  \item{pThreshold}{ threshold of Fisher's exact test for considering modules to have a significant
overlap. }
}
\details{
  Each column of \code{source} is treated separately. Source and reference labels are assumed to be of
the same type, that is noth should be either numeric or character. If the labels are character, they are
assumed to be color labels such as the ones returned by \code{\link{standardColors}}. 

  The function calculates the overlap of the \code{source} and \code{reference} modules using Fisher's
exact test. It then attempts to relabel \code{source} modules such that modules with the highest overlap
with the \code{reference} modules have the same color. Where this is not possible (for example because
one reference module has the highest overlap with two source modules), the source modules will be
relabeled using labels that are not present among the reference labels.
}
\value{
  A vector (if the input \code{source} labels are a vector) or a matrix (if the input \code{source}
labels are a matrix) of the new labels.
}

\author{ Peter Langfelder }
\seealso{ \code{\link{standardColors}} }
\keyword{ misc }
