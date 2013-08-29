\name{exportNetworkToCytoscape}
\alias{exportNetworkToCytoscape}
\title{ Export network to Cytoscape }
\description{
  This function exports a network in edge and node list files in a format suitable for importing to
Cytoscape.
}
\usage{
exportNetworkToCytoscape(
   adjMat,
   edgeFile = NULL,
   nodeFile = NULL,
   weighted = TRUE,
   threshold = 0.5,
   nodeNames = NULL,
   altNodeNames = NULL,
   nodeAttr = NULL,
   includeColNames = TRUE)
}
\arguments{
  \item{adjMat}{ adjacency matrix giving connection strengths among the nodes in the network. }
  \item{edgeFile}{ file name of the file to contain the edge information. }
  \item{nodeFile}{ file name of the file to contain the node information. }
  \item{weighted}{ logical: should the exported network be weighted? }
  \item{threshold}{ adjacency threshold for including edges in the output.  }
  \item{nodeNames}{ names of the nodes. If not given, \code{dimnames} of \code{adjMat} will be used. }
  \item{altNodeNames}{ optional alternate names for the nodes, for example gene names if nodes are
labeled by probe IDs. } 
  \item{nodeAttr}{ optional node attribute, for example module color. Can be a vector or a data frame. }
  \item{includeColNames}{ logical: should column names be included in the output files? Note that
Cytoscape can read files both with and without column names.  }
}
\details{
  If the corresponding file names are supplied, the edge and node data is written to the appropriate
files. The edge and node data is also returned as return value (see below). }
\value{
  A list with the following componens:
  \item{egdeData}{a data frame containing the edge data, with one row per edge}
  \item{nodeData}{a data frame containing the node data, with one row per node}
}
\author{ Peter Langfelder}
\seealso{ \code{\link{exportNetworkToVisANT}}}

\keyword{ misc }
