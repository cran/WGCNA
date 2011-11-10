\name{labels2colors}
\alias{labels2colors}
\title{Convert numerical labels to colors. }
\description{
  Converts a vector or array of numerical labels into a corresponding vector or array of colors corresponding to the labels. 
}
\usage{
labels2colors(labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey",
              commonColorCode = TRUE)
}
\arguments{
  \item{labels}{Vector or matrix of non-negative integer or other (such as character) labels. See details.}
  \item{zeroIsGrey}{If TRUE, labels 0 will be assigned color grey. Otherwise, labels below 1 will trigger
an error.}
  \item{colorSeq}{Color sequence corresponding to labels. If not given, a standard sequence will be
used.}
  \item{naColor}{Color that will encode missing values. }
  \item{commonColorCode}{logical: if \code{labels} is a matrix, should each column have its own colors? }
}
\details{
If \code{labels} is numeric, it is used directly as index to the standard color sequence. If 0 is present
among the labels and \code{zeroIsGrey=TRUE}, labels 0 are given grey color.

If \code{labels} is not numeric, its columns are turned into factors and the numeric representation of each
factor is used to assign the corresponding colors. In this case \code{commonColorCode} governs whether each
column gets its own color code, or whether the color code will be universal.

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
labels = matrix(letters[1:9], 3,3);
labels2colors(labels)
# Note the difference when commonColorCode = FALSE
labels2colors(labels, commonColorCode = FALSE)
}

\keyword{color}
