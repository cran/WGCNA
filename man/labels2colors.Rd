\name{labels2colors}
\alias{labels2colors}
\title{Convert numerical labels to colors. }
\description{
  Converts a vector or array of numerical labels into a corresponding vector or array of colors corresponding to the labels. 
}
\usage{
labels2colors(labels, zeroIsGrey = TRUE, colorSeq = NULL)
}
\arguments{
  \item{labels}{Vector of non-negative integer labels.}
  \item{zeroIsGrey}{If TRUE, labels 0 will be assigned color grey. Otherwise, labels below 1 will trigger
an error.}
  \item{colorSeq}{Color sequence corresponding to labels. If not given, a standard sequence will be
used.}
}
\details{
The standard sequence start with well-distinguishable colors, and after about 40 turns into a
quasi-random sampling of all colors available in R with the exception of all shades of grey (and gray).

If the input \code{labels} have a dimension attribute, it is copied into the output, meaning the
dimensions of the returned value are the same as those of the input \code{labels}.
}
\value{
  A vector or array of character strings of the same length or dimensions as \code{labels}.
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\examples{
labels = c(0:20);
labels2colors(labels);
}

\keyword{color}
