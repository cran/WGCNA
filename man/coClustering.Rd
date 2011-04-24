\name{coClustering}
\alias{coClustering}
\title{
  Co-clustering measure of cluster preservation between two clusterings
}
\description{
The function calculates the co-clustering statistics for each module in the reference clustering.
}
\usage{
coClustering(clusters.ref, clusters.test, tupletSize = 2, unassignedLabel = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters.ref}{
  Reference input clustering. A vector in which each element gives the cluster label of an object.
}
  \item{clusters.test}{
  Test input clustering. Must be a vector of the same size as \code{cluster.ref}. 
}
  \item{tupletSize}{
  Co-clutering tuplet size.
}
  \item{unassignedLabel}{
  Optional specification of a clustering label that denotes unassigned objects. Objects with this label are
excluded from the calculation.
}
}
\details{
Co-clustering of cluster q in the reference clustering and cluster q' in the test clustering measures the
overlap of clusters q and q' by the number of tuplets that can be chosen from the overlap of clusters q and
q' relative to the number of tuplets in cluster q. To arrive at a co-clustering measure for cluster q, we
sum the co-clustering of q and q' over all clusters q' in the test clustering. A value close to 1 indicates
high preservation of the reference cluster in the test clustering, while a value close to zero indicates a
low preservation.
}
\value{
A vector in which each component corresponds to a cluster in the reference clustering. Entries give the
co-clustering measure of cluster preservation. 
}
\references{
For example, see Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is My Network Module Preserved and
Reproducible? PLoS Comput Biol 7(1): e1001057. Co-clustering is discussed in the Methods Supplement
(Supplementary text 1) of that article. 
}
\author{
 Peter Langfelder
}
\seealso{
\code{\link{modulePreservation}} for a large suite of module preservation statistics
\code{\link{coClustering.permutationTest}} for a permutation test for co-clustering significance
}
\examples{

  # An example with random (unrelated) clusters:

  set.seed(1);
  nModules = 10;
  nGenes = 1000;
  cl1 = sample(c(1:nModules), nGenes, replace = TRUE);
  cl2 = sample(c(1:nModules), nGenes, replace = TRUE);
  coClustering(cl1, cl2)

  # For the same reference and test clustering:

  coClustering(cl1, cl1)


}
\keyword{misc}
