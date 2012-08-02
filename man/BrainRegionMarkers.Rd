\name{BrainRegionMarkers}
\alias{BrainRegionMarkers}
\docType{data}
\title{Gene Markers for Regions of the Human Brain}
\description{
This matrix gives a predefined set of marker genes for many regions of the human brain, using data from the Allen Human Brain Atlas (http://human.brain-map.org/) as reported in: Hawrylycz MJ, Lein ES, Guillozet-Bongaarts AL, Shen EH, Ng L, Miller JA, et al. (2012) An Anatomically Comprehensive Atlas of the Adult Human Brain Transcriptome. Nature (in press).  It is used with userListEnrichment to search user-defined gene lists for enrichment.
}
\usage{data(BrainRegionMarkers)}
\format{
A 28477 x 2 matrix of characters containing Gene / Category pairs.  The first column (Gene) lists genes corresponding to a given category (second column).  Each Category entry is of the form <Brain Region>_<Marker Type>__HBA.  Note that the matrix is sorted first by Category and then by Gene, such that all genes related to the same category are listed sequentially.
}
\source{
For references used in this variable, or other information, please see \code{\link{userListEnrichment}}
}
\examples{
data(BrainRegionMarkers)
head(BrainRegionMarkers)
}
\keyword{datasets}
