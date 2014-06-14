\name{mtd.simplify}
\alias{mtd.simplify}

\title{
If possible, simplify a multiData structure to a 3-dimensional array.
}
\description{
This function attempts to put all \code{data} components into a 3-dimensional array, with the last
dimension corresponding to the sets. This is only possible if all \code{data} components are matrices or
data frames with the same dimensiosn. 
}
\usage{
mtd.simplify(multiData)
}
\arguments{
  \item{multiData}{
A multiData structure in the "strict" sense (see below).
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

This function assumes a "strict" multiData structure.
}
\value{
A 3-dimensional array collecting all \code{data} components.
}
\author{
Peter Langfelder
}
\note{
The function is relatively fragile and may fail. Use at your own risk.
}

\seealso{
\code{\link{multiData}} to create a multiData structure;

\code{\link{multiData2list}} for converting multiData structures to plain lists.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
