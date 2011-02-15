\name{plotMEpairs}
\alias{plotMEpairs}
\title{ Pairwise scatterplots of eigengenes}
\description{
  The function produces a matrix of plots containing pairwise scatterplots of given eigengenes, the
distribution of their values and their pairwise correlations.
}
\usage{
plotMEpairs(
   datME, 
   y = NULL, 
   main = "Relationship between module eigengenes", 
   clusterMEs = TRUE, 
   ...)
}
\arguments{
  \item{datME}{  a data frame containing expression data, with rows corresponding to samples and columns
to genes. Missing values are allowed and will be ignored. }

  \item{y}{ optional microarray sample trait vector. Will be treated as an additional eigengene. }
  \item{main}{ main title for the plot. }
  \item{clusterMEs}{ logical: should the module eigengenes be ordered by their dendrogram? }
  \item{\dots}{ additional graphical parameters to the function \code{\link{pairs}} }
}
\details{
  The function produces an NxN matrix of plots, where N is the number of eigengenes. In the upper
traingle it plots pairwise scatterplots of module eigengenes (plus the trait \code{y}, if given). On the
diagonal it plots histograms of sample values for each eigengene. Below the diagonal, it displays the
pairwise correlations of the eigengenes.
}
\value{
  None.
}
\author{ Steve Horvath }
\seealso{ \code{\link{pairs}} }
\keyword{ hplot }% __ONLY ONE__ keyword per line
