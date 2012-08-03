\name{BloodLists}
\alias{BloodLists}
\docType{data}
\title{Blood Cell Types with Corresponding Gene Markers}
\description{
This matrix gives a predefined set of marker genes for many blood cell types, as reported in several previously-published studies.  It is used with userListEnrichment to search user-defined gene lists for enrichment.
}
\usage{data(BloodLists)}
\format{
A 2048 x 2 matrix of characters containing Gene / Category pairs.  The first column (Gene) lists genes corresponding to a given category (second column).  Each Category entry is of the form <Blood cell type>__<reference>, where the references can be found at \code{\link{userListEnrichment}}.  Note that the matrix is sorted first by Category and then by Gene, such that all genes related to the same category are listed sequentially.
}
\source{
For references used in this variable, please see \code{\link{userListEnrichment}}
}
\examples{
data(BloodLists)
head(BloodLists)
}
\keyword{datasets}
