\name{sizeRestrictedClusterMerge}
\alias{sizeRestrictedClusterMerge}
\title{
Cluter merging with size restrictions
}
\description{
This function merges clusters by correlation of the first principal components such that the resulting merged clusters do
not exceed a given maximum size.
}
\usage{
sizeRestrictedClusterMerge(
    datExpr,
    clusters,
    clusterSizes = NULL,
    centers = NULL,
    maxSize,
    networkType = "unsigned",
    verbose = 0,
    indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{
Data on which the clustering is based (e.g., expression data). Variables are in columns and observations (samples) in
rows.
}
  \item{clusters}{
A vector with element per variable (column) in \code{datExpr} giving the cluster label for the corresponding variable.
}
  \item{clusterSizes}{
Optional pre-calculated cluster sizes. If not given, will be determined from given \code{clusters}.
}
  \item{centers}{
Optional pre-calculaed cluster centers (first principal components/singular vectors). 
If not given, will be calculated from given data and
cluster assignments.
}
  \item{maxSize}{
Maximum allowed size of merged clusters. If any of the given \code{clusters} are larger than \code{maxSize}, they will not
be changed.
}
  \item{networkType}{
One of \code{"unsigned"} and \code{"signed"}. Determines whether clusters with negatively correlated representatives will
be considered similar (\code{"unsigned"}) or dis-similar (\code{"signed"}). 
}
  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
The function iteratively merges two closest clusters subject to the constraint that the merged cluster size cannot exceed
maxSize. Merging stops when no two clusters can be merged without exceeding the maximum size.
}
\value{
A list with two components
\item{clusters}{A numeric vector with one component per input gene, giving the cluster number in
which the gene is assigned. }
\item{centers}{Cluster centers, that is their first principal components/singular vectors. }
}
\author{
Peter Langfelder
}
\seealso{
The last step in \code{\link{projectiveKMeans}} uses this function.
}
\keyword{cluster}
