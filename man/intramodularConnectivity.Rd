\name{intramodularConnectivity}
\alias{intramodularConnectivity}
\title{ Calculation of intramodular connectivity }
\description{
  Calculates intramodular connectivity, i.e., connectivity of nodes to other nodes within the same
module.
}
\usage{
intramodularConnectivity(adjMat, colors, scaleByMax = FALSE)
}
\arguments{
  \item{adjMat}{ adjacency matrix, a square, symmetric matrix with entries between 0 and 1. }
  \item{colors}{ module labels. A vector of length \code{ncol(adjMat)} giving a module label for each
gene (node) of the network.}
  \item{scaleByMax}{ logical: should intramodular connectivities be scaled by the maximum IM connectivity
in each module? }
}
\details{
  The module labels can be numeric or character. For each node (gene), the function sums adjacency
entries (excluding the diagonal) to other nodes within the same module. Optionally, the connectivities
can be scaled by the maximum connectivy in each module. 
}
\value{
  A data frame with 4 columns giving the total connectivity, intramodular connectivity, extra-modular
connectivity, and the difference of the intra- and extra-modular connectivities for all genes.
}
\references{ Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology
2007, 1:24 }

\author{ Steve Horvath and Peter Langfelder}

\seealso{ \code{\link{adjacency}} }

\keyword{ misc }% __ONLY ONE__ keyword per line
