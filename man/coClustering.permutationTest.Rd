\name{coClustering.permutationTest}
\alias{coClustering.permutationTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Permutation test for co-clustering
}
\description{
 This function calculates permutation Z statistics that measure how different the co-clustering of modules
in a reference and test clusterings is from random. 
}
\usage{
coClustering.permutationTest(
      clusters.ref, clusters.test, 
      tupletSize = 2, 
      nPermutations = 100, 
      unassignedLabel = 0, 
      randomSeed = 12345, verbose = 0, indent = 0)
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
  \item{nPermutations}{
   Number of permutations to execute. Since the function calculates parametric p-values, a relatively small
number of permutations (at least 50) should be sufficient.
}
  \item{unassignedLabel}{
  Optional specification of a clustering label that denotes unassigned objects. Objects with this label are
excluded from the calculation.
}
  \item{randomSeed}{
  Random seed for initializing the random number generator. If \code{NULL}, the generator is not initialized
(useful for calling the function sequentially). The default assures reproducibility. 
}
  \item{verbose}{
  If non-zero, function will print out progress messages.
}
  \item{indent}{
  Indentation for progress messages. Each unit adds two spaces.
}
}
\details{
This function performs a permutation test to determine whether observed co-clustering statistics are
significantly different from those expected by chance. It returns the observed co-clustering as well as the
permutation Z statistic, calculated as \code{(observed - mean)/sd}, where \code{mean} and \code{sd} are the
mean and standard deviation of the co-clustering when the test clustering is repeatedly randomly permuted. 
}
\value{
  \item{observed }{the observed co-clustering measures for clusters in \code{clusters.ref} }
  \item{Z}{permutation Z statics}
  \item{permuted.mean}{means of the co-clustering measures when the test clustering is permuted}
  \item{permuted.sd}{standard deviations of the co-clustering measures when the test clustering is permuted}
  \item{permuted.cc}{values of the co-clustering measure for each permutation of the test clustering. A
matrix of dimensions (number of permutations)x(number of clusters in reference clustering). }

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
\code{\link{coClustering}} for calculation of the "observed" co-clustering measure 
\code{\link{modulePreservation}} for a large suite of module preservation statistics
}

\examples{
  set.seed(1);
  nModules = 5;
  nGenes = 100;
  cl1 = sample(c(1:nModules), nGenes, replace = TRUE);
  cl2 = sample(c(1:nModules), nGenes, replace = TRUE);
  
  cc = coClustering(cl1, cl2)

  # Choose a low number of permutations to make the example fast
  ccPerm = coClustering.permutationTest(cl1, cl2, nPermutations = 20, verbose = 1);

  ccPerm$observed
  ccPerm$Z

  # Combine cl1 and cl2 to obtain clustering that is somewhat similar to cl1:

  cl3 = cl2;
  from1 = sample(c(TRUE, FALSE), nGenes, replace = TRUE);
  cl3[from1] = cl1[from1];

  ccPerm = coClustering.permutationTest(cl1, cl3, nPermutations = 20, verbose = 1);

  # observed co-clustering is higher than before:
  ccPerm$observed

  # Note the high preservation Z statistics:
  ccPerm$Z
}
\keyword{misc}
