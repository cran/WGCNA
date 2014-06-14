\name{mtd.setAttr}
\alias{mtd.setAttr}

\title{
Set attributes on each component of a multiData structure
}
\description{
Set attributes on each \code{data} component of a multiData structure
}
\usage{
mtd.setAttr(multiData, attribute, valueList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiData}{
A multiData structure.
}
  \item{attribute}{
Name for the attribute to be set
}
  \item{valueList}{
List that gives the attribute value for each set in the multiData structure.
}
}
\value{
The input multiData with the attribute set on each \code{data} component.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{multiData}} to create a multiData structure;

\code{isMultiData} for a description of the multiData structure.
}
\keyword{misc}% __ONLY ONE__ keyword per line
