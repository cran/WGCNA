\name{ImmunePathwayLists}
\alias{ImmunePathwayLists}
\docType{data}
\title{Immune Pathways with Corresponding Gene Markers}
\description{
This matrix gives a predefined set of marker genes for many immune response pathways, as assembled by Brian Modena (a member of Daniel R Salomon's lab at Scripps Research Institute), and colleagues.  It is used with userListEnrichment to search user-defined gene lists for enrichment.
}
\usage{data(ImmunePathwayLists)}
\format{
A 3597 x 2 matrix of characters containing Gene / Category pairs.  The first column (Gene) lists genes corresponding to a given category (second column).  Each Category entry is of the form <Immune Pathway>__ImmunePathway.  Note that the matrix is sorted first by Category and then by Gene, such that all genes related to the same category are listed sequentially.
}
\source{
For more information about this list, please see \code{\link{userListEnrichment}}
}
\examples{
data(ImmunePathwayLists)
head(ImmunePathwayLists)
}
\keyword{datasets}
