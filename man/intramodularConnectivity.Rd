\name{intramodularConnectivity}
\alias{intramodularConnectivity}
\alias{intramodularConnectivity.fromExpr}
\title{ Calculation of intramodular connectivity }
\description{
  Calculates intramodular connectivity, i.e., connectivity of nodes to other nodes within the same
module.
}
\usage{
intramodularConnectivity(adjMat, colors, scaleByMax = FALSE)

intramodularConnectivity.fromExpr(datExpr, colors, 
              corFnc = "cor", corOptions = "use = 'p'",
              weights = NULL,
              distFnc = "dist", distOptions = "method = 'euclidean'",
              networkType = "unsigned", power = if (networkType=="distance") 1 else 6,
              scaleByMax = FALSE,
              ignoreColors = if (is.numeric(colors)) 0 else "grey",
              getWholeNetworkConnectivity = TRUE)

}
\arguments{
  \item{adjMat}{ adjacency matrix, a square, symmetric matrix with entries between 0 and 1. }
  \item{colors}{ module labels. A vector of length \code{ncol(adjMat)} giving a module label for each
gene (node) of the network.}
  \item{scaleByMax}{ logical: should intramodular connectivities be scaled by the maximum IM connectivity
in each module? }

  \item{datExpr}{ data frame or matrix containing expression data. Columns correspond to genes and rows to
samples.}

  \item{corFnc}{ character string specifying the function to be used to calculate co-expression
similarity for correlation networks.
Defaults to Pearson correlation. Any function returning values between -1 and 1 can be used. }

  \item{corOptions}{ character string specifying additional arguments to be passed to the function given
by \code{corFnc}. Use \code{"use = 'p', method = 'spearman'"} to obtain Spearman correlation.   }

  \item{weights}{optional matrix of the same dimensions as \code{datExpr}, giving the weights for individual
observations in \code{datExpr}. These will be passed on to the correlation function.}

  \item{distFnc}{ character string specifying the function to be used to calculate co-expression
similarity for distance networks. Defaults to the function \code{\link{dist}}.
Any function returning non-negative values can be used.}

  \item{distOptions}{ character string specifying additional arguments to be passed to the function given
by \code{distFnc}. For example, when the function  \code{\link{dist}} is used, the argument \code{method}
can be used to specify various ways of computing the distance. }

  \item{networkType}{network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}, \code{"distance"}. }

  \item{power}{soft thresholding power. }

  \item{ignoreColors}{level(s) of \code{colors} that identifies unassigned genes. The intramodular
connectivity in this "module" will not be calculated.} 

  \item{getWholeNetworkConnectivity}{logical: should whole-network connectivity be computed as well? For
large networks, this can be quite time-consuming.}



}
\details{
  The module labels can be numeric or character. For each node (gene), the function sums adjacency
entries (excluding the diagonal) to other nodes within the same module. Optionally, the connectivities
can be scaled by the maximum connectivy in each module. 
}
\value{
  If input \code{getWholeNetworkConnectivity} is \code{TRUE}, a data frame with 4 columns giving the total connectivity, intramodular connectivity, extra-modular
connectivity, and the difference of the intra- and extra-modular connectivities for all genes; otherwise a
vector of intramodular connectivities, 
}
\references{ Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology
2007, 1:24 }

\author{ Steve Horvath and Peter Langfelder}

\seealso{ \code{\link{adjacency}} }

\keyword{ misc }% __ONLY ONE__ keyword per line
