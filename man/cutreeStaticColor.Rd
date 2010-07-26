\name{cutreeStaticColor}
\alias{cutreeStaticColor}
\title{ Constant height tree cut using color labels}
\description{
  Cluster detection by a constant height cut of a hierarchical clustering dendrogram.
}
\usage{
cutreeStaticColor(dendro, cutHeight = 0.9, minSize = 50)
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
  A character vector giving color labels of objects, with "grey" meaning unassigned. The largest cluster
is conventionally labeled "turquoise", next "blue" etc. Run \code{standardColors()} to see the sequence
of standard color labels.
}
\author{ Peter Langfelder }
\seealso{ \code{\link{hclust}} for hierarchical clustering, \code{\link{cutree}} and
\code{\link{cutreeStatic}} for other constant-height branch cuts, \code{\link{standardColors}} to see
the sequence of color labels that can be assigned.}
\keyword{misc}
