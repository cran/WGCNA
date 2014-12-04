\name{branchSplitFromStabilityLabels}
\alias{branchSplitFromStabilityLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Branch split (dissimilarity) statistic derived from labels determined from a stability study
}
\description{
This function evaluates how different two branches are based on a series of cluster labels that are usually
obtained in a stability study but can in principle be arbitrary. The idea is to quantify how well
membership on the two tested branches can be predicted from clusters in the given stability labels.
}
\usage{
branchSplitFromStabilityLabels(
   branch1, branch2, 
   stabilityLabels, 
   ignoreLabels = 0,
   ...)
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
indistinguishable (based on the test clustering). 

Formally, for each cluster C in each clustering in \code{stabilityLabels}, 
its contribution to the branch similarity
is min(r1, r2), where r1 = |intersect(C, branch1)|/|branch1| and r2 = |intersect(C, branch2)|/|branch2|.
The statistics for clusters in each clustering are added; the sums are then averaged across the
clusterings. Since the result is a similarity statistic, the final dissimilarity is defined as
1-similarity. The dissimilarity ranges between 0 (branch1 and branch2 are indistinguishable) and 1 (branch1
and branch2 are perfectly distinguishable).

This is a very simple statistic that does not attempt to correct for the similarity that would be expected
by chance. 
}
\value{
Branch dissimilarity (a single number between 0 and 1).
}
\author{
Peter Langfelder
}

\seealso{
This function is utilized in \code{\link{blockwiseModules}} and \code{\link{blockwiseConsensusModules}}.
}
\keyword{misc}
