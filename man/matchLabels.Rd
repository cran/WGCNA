\name{matchLabels}
\alias{matchLabels}
\title{ Relabel module labels to best match the given reference labels }
\description{
  Given a \code{source} and  \code{reference} vectors of module labels, the function produces a module
labeling that is equivalent to \code{source}, but individual modules are re-labeled so that modules with
significant overlap in \code{source} and \code{reference} have the same labels.
}
\usage{
matchLabels(source, 
            reference, 
            pThreshold = 5e-2,
            na.rm = TRUE,
            ignoreLabels = if (is.numeric(reference)) 0 else "grey",
            extraLabels = if (is.numeric(reference)) c(1:1000) else standardColors()
            )
}
\arguments{
  \item{source}{ a vector or a matrix of reference labels. The labels may be numeric or character. }
  \item{reference}{ a vector of reference labels. }
  \item{pThreshold}{ threshold of Fisher's exact test for considering modules to have a significant
overlap. }
  \item{na.rm}{logical: should missing values in either \code{source} or \code{reference} be removed? If
not, missing values may be treated as a standard label or the function may throw an error 
(exact behaviour depends on whether the input labels
are numeric or not).}
  \item{ignoreLabels}{labels in \code{source} and \code{reference} to be considered unmatchable. These
labels are excluded from the re-labeling procedure.}
  \item{extraLabels}{a vector of labels for modules in \code{source} that cannot be matched to any modules
in \code{reference}. The user should ensure that this vector contains enough labels since the function
automatically removes a values that occur in either \code{source}, \code{reference} or \code{ignoreLabels},
to avoid possible confusion. }
}
\details{
  Each column of \code{source} is treated separately. Unlike in previous version of this function, source
and reference labels can be any labels, not necessarily of the same type.

  The function calculates the overlap of the \code{source} and \code{reference} modules using Fisher's
exact test. It then attempts to relabel \code{source} modules such that each \code{source} module gets the
label of the \code{reference} module that it overlaps most with, subject to not renaming two \code{source}
modules to the same \code{reference} module. (If two \code{source} modules point to the same
\code{reference} module, the one with the more significant overlap is chosen.)

Those \code{source} modules that cannot be matched to a \code{reference} module are labeled using
those labels from \code{extraLabels} that do not occur in either of \code{source}, \code{reference} or
\code{ignoreLabels}.

}
\value{
  A vector (if the input \code{source} labels are a vector) or a matrix (if the input \code{source}
labels are a matrix) of the new labels.
}

\author{ Peter Langfelder }
\seealso{ 

\code{\link{overlapTable}} for calculation of overlap counts and p-values;

\code{\link{standardColors}} for standard non-numeric WGCNA labels.
}
\keyword{ misc }
