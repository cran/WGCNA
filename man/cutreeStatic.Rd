\name{cutreeStatic}
\alias{cutreeStatic}
\title{ Constant-height tree cut } 
\description{
  Module detection in hierarchical dendrograms using a constant-height tree cut. Only branches whose size
is at least \code{minSize} are retained. 
}
\usage{
cutreeStatic(dendro, cutHeight = 0.9, minSize = 50)
}
\arguments{
  \item{dendro}{ a hierarchical clustering dendrogram such as returned by \code{\link{hclust}}. }
  \item{cutHeight}{ height at which branches are to be cut. }
  \item{minSize}{ minimum number of object on a branch to be considered a cluster. }
}
\details{
  This function performs a straightforward constant-height cut as implemented by \code{\link{cutree}},
then calculates the number of objects on each branch and only keeps branches that have at least
\code{minSize} objects on them.
}
\value{
  A numeric vector giving labels of objects, with 0 meaning unassigned. The largest cluster
is conventionally labeled 1, the next largest 2, etc. 
}
\author{ Peter Langfelder }
\seealso{ \code{\link{hclust}} for hierarchical clustering, \code{\link{cutree}} and
\code{\link{cutreeStatic}} for other constant-height branch cuts, \code{\link{standardColors}} to convert
the retuned numerical lables into colors for easier visualization. }
\keyword{misc}

