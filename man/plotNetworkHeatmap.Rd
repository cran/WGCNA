\name{plotNetworkHeatmap}
\alias{plotNetworkHeatmap}
\title{ Network heatmap plot }
\description{
  Network heatmap plot.
}
\usage{
plotNetworkHeatmap(
  datExpr, 
  plotGenes, 
  weights = NULL,
  useTOM = TRUE, 
  power = 6, 
  networkType = "unsigned", 
  main = "Heatmap of the network")
}
\arguments{
  \item{datExpr}{ a data frame containing expression data, with rows corresponding to samples and columns
to genes. Missing values are allowed and will be ignored. }

  \item{plotGenes}{ a character vector giving the names of genes to be included in the plot. The names
will be matched against \code{names(datExpr)}. }

  \item{weights}{optional observation weights for \code{datExpr} to be used in correlation calculation.
A matrix of the same dimensions as \code{datExpr}, containing non-negative weights. Only used with Pearson
correlation.}

  \item{useTOM}{ logical: should TOM be plotted (\code{TRUE}), or correlation-based adjacency
(\code{FALSE})? }

  \item{power}{ soft-thresholding power for network construction. }

  \item{networkType}{ a character string giving the newtork type. Recognized values are (unique
abbreviations of) \code{"unsigned"}, \code{"signed"}, and \code{"signed hybrid"}. }

  \item{main}{ main title for the plot. }
}
\details{
  The function constructs a network from the given expression data (selected by \code{plotGenes}) using
the soft-thresholding procedure, optionally calculates Topological Overlap (TOM) and plots a heatmap of
the network. 

Note that all network calculations are done in one block and may fail due to memory allocation issues for
large numbers of genes. 

}
\value{
  None.
}
\references{ 
Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
}
\author{ Steve Horvath }
\seealso{ \code{\link{adjacency}}, \code{\link{TOMsimilarity}} }
\keyword{ hplot }% __ONLY ONE__ keyword per line
