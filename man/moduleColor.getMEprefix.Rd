\name{moduleColor.getMEprefix}
\alias{moduleColor.getMEprefix}
\title{Get the prefix used to label module eigengenes.}
\description{
  Returns the currently used prefix used to label module eigengenes. 
When returning module eigengenes in a dataframe,
names of the corresponding columns will start with the given prefix.
}
\usage{
moduleColor.getMEprefix()
}
\details{
  Returns the prefix used to label module eigengenes. When returning module eigengenes in a dataframe,
names of the corresponding columns will consist of the corresponfing color label preceded by the
given prefix. For example, if the prefix is "PC" and the module is turquoise, the corresponding
module eigengene will be labeled "PCturquoise". Most of old code assumes "PC", but "ME" is more
instructive and used in some newer analyses.
}
\value{
  A character string.
}
\note{
  Currently the standard prefix is \code{"ME"} and there is no way to change it.
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{
 \code{\link{moduleEigengenes}}
}
\keyword{misc}
