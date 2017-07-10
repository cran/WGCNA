\name{branchSplitFromStabilityLabels}
\alias{branchSplitFromStabilityLabels}
\alias{branchSplitFromStabilityLabels.individualFraction}
\alias{branchSplitFromStabilityLabels.prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Branch split (dissimilarity) statistics derived from labels determined from a stability study
}
\description{
These functions evaluate how different two branches are based on a series of cluster labels that are usually
obtained in a stability study but can in principle be arbitrary. The idea is to quantify how well
membership on the two tested branches can be predicted from clusters in the given stability labels.
}
\usage{
branchSplitFromStabilityLabels(
   branch1, branch2, 
   stabilityLabels, 
   ignoreLabels = 0,
   ...)

branchSplitFromStabilityLabels.prediction(
            branch1, branch2,
            stabilityLabels, ignoreLabels = 0, ...)

branchSplitFromStabilityLabels.individualFraction(
            branch1, branch2,
            stabilityLabels, ignoreLabels = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{branch1}{
A vector of indices giving members of branch 1.
}
  \item{branch2}{
A vector of indices giving members of branch 1.
}
  \item{stabilityLabels}{
A matrix of cluster labels. Each column corresponds to one clustering and each row to one object (whose
indices \code{branch1} and \code{branch2} refer to). 
}
  \item{ignoreLabels}{
Label or labels that do not constitute proper clusters in \code{stabilityLabels}, for example because they
label unassigned objects. 
}
  \item{\dots}{
Ignored.
}
}
\details{
The idea is to measure how well clusters in \code{stabilityLabels} can distinguish the two given branches.
For example, if a cluster C intersects with branch1 but not branch2, it can distinguish branches 1 and 2
perfectly. On the other hand, if there is a cluster C that contains both branch 1 and branch 2, 
the two branches are
indistinguishable (based on the test clustering).  The three functions differ in the details of the
similarity calculation.  

\code{branchSplitFromStabilityLabels.individualFraction}: Currently the recommended branch split calculation
method, and default for \code{\link{hierarchicalConsensusModules}}. 
For each branch and all clusters that overlap with the branch (not necessarily with the other
branch), calculate the fraction of the cluster objects (restricted to the two
branches) that belongs to the branch. For each branch, sum these fractions over all clusters. 
If this number is relatively low, around 0.5, it means most elements are in
non-discriminative clusters.

\code{branchSplitFromStabilityLabels}: This was the original branch split measure and for backward
compatibility it still is the default method in \code{\link{blockwiseModules}} and
\code{\link{blockwiseConsensusModules}}. For each cluster C in each clustering in \code{stabilityLabels}, 
its contribution to the branch similarity
is min(r1, r2), where r1 = |intersect(C, branch1)|/|branch1| and r2 = |intersect(C, branch2)|/|branch2|.
The statistics for clusters in each clustering are added; the sums are then averaged across the
clusterings. 

\code{branchSplitFromStabilityLabels.prediction}: Use only for experiments, not recommended for actual
analyses because it is not stable under small changes in the branch membership.
For each cluster that overlaps with both branches, count the objects in the branch with which the cluster 
has a smaller overlap
and add it to the score for that branch. The final counts divided by number of genes on branch give a
"indistinctness" score; take the larger of the two indistinctness scores and call this the similarity.

Since the result of the last two calculations is a similarity statistic, the final dissimilarity is defined as
1-similarity. The dissimilarity ranges between 0 (branch1 and branch2 are indistinguishable) and 1 (branch1
and branch2 are perfectly distinguishable).

These statistics are quite simple and do not correct for similarity that would be expected by chance. On the
other hand, all 3 statistics are fairly (though not perfectly) stable under splitting and joining of clusters
in \code{stabilityLabels}.
}
\value{
Branch dissimilarity (a single number between 0 and 1).
}
\author{
Peter Langfelder
}

\seealso{
These function are utilized in \code{\link{blockwiseModules}}, \code{\link{blockwiseConsensusModules}} and
\code{\link{hierarchicalConsensusModules}}.
}
\keyword{misc}
