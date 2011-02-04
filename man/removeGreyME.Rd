\name{removeGreyME}
\alias{removeGreyME}
\title{Removes the grey eigengene from a given collection of eigengenes. }
\description{
  Given module eigengenes either in a single data frame or in a multi-set format, removes the grey
eigengenes from each set. If the grey eigengenes are not found, a warning is issued.
}
\usage{
removeGreyME(MEs, greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
}
\arguments{
  \item{MEs}{Module eigengenes, either in a single data frame (typicaly for a single set), or in a
multi-set format. See \code{\link{checkSets}} for a description of the multi-set format.}
  \item{greyMEName}{Name of the module eigengene (in each corresponding data frame) that corresponds to
the grey color. This will typically be "PCgrey" or "MEgrey". If the module eigengenes were calculated
using standard functions in this library, the default should work.}
}
\value{
  Module eigengenes in the same format as input (either a single data frame or a vector of lists) with the
grey eigengene removed. 
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}

\keyword{misc}
