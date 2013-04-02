\name{PWLists}
\alias{PWLists}
\docType{data}
\title{Pathways with Corresponding Gene Markers - Compiled by Mike Palazzolo and Jim Wang from CHDI}
\description{
This matrix gives a predefined set of marker genes for many immune response pathways, as assembled by Mike Palazzolo and Jim Wang from CHDI, and colleagues.  It is used with userListEnrichment to search user-defined gene lists for enrichment.
}
\usage{data(PWLists)}
\format{
A 124350 x 2 matrix of characters containing 2724 Gene / Category pairs.  The first column (Gene) lists genes corresponding to a given category (second column).  Each Category entry is of the form <gene set>__<reference>.  
}
\source{
For more information about this list, please see \code{\link{userListEnrichment}}
}
\examples{
data(PWLists)
head(PWLists)
}
\keyword{datasets}
