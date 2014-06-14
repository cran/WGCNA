\name{mtd.setColnames}
\alias{mtd.setColnames}
\alias{mtd.colnames}
\title{
Get and set column names in a multiData structure.
}
\description{
Get and set column names on each \code{data} component in a multiData structure.

}
\usage{
mtd.colnames(multiData)
mtd.setColnames(multiData, colnames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiData}{
A multiData structure
}
  \item{colnames}{
A vector (coercible to character) of column names. 
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

The \code{mtd.colnames} and \code{mtd.setColnames} assume (and checks for) a "strict" multiData structure.

}
\value{
\code{mtd.colnames} returns the vector of column names of the \code{data} component. The function assumes
the column names in all sets are the same. 

\code{mtd.setColnames} returns the multiData structure with the column names set in all \code{data}
components.
 }
\author{
Peter Langfelder
}

\seealso{
\code{\link{multiData}} to create a multiData structure.
}

\keyword{ misc }% __ONLY ONE__ keyword per line
