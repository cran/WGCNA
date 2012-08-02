\name{SCsLists}
\alias{SCsLists}
\docType{data}
\title{Stem Cell-Related Genes with Corresponding Gene Markers}
\description{
This matrix gives a predefined set of genes related to several stem cell (SC) types, as reported in two previously-published studies.  It is used with userListEnrichment to search user-defined gene lists for enrichment.
}
\usage{data(SCsLists)}
\format{
A 14003 x 2 matrix of characters containing Gene / Category pairs.  The first column (Gene) lists genes corresponding to a given category (second column).  Each Category entry is of the form <Stem cell-related category>__<reference>, where the references can be found at \code{\link{userListEnrichment}}.  Note that the matrix is sorted first by Category and then by Gene, such that all genes related to the same category are listed sequentially.
}
\source{
For references used in this variable, please see \code{\link{userListEnrichment}}
}
\examples{
data(SCsLists)
head(SCsLists)
}
\keyword{datasets}
