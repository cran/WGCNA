\name{fundamentalNetworkConcepts}
\alias{fundamentalNetworkConcepts}
\title{ Calculation of fundamental network concepts from an adjacency matrix. }
\description{
This function computes fundamental network concepts (also known as network indices or statistics) based
on an adjacency matrix and optionally a node significance measure. These network concepts are defined for
any symmetric adjacency matrix (weighted and unweighted). The network concepts are described in Dong and
Horvath (2007) and Horvath and Dong (2008). 
Fundamental network concepts are defined as a function of the off-diagonal elements of an adjacency
matrix \code{adj} and/or a node significance measure \code{GS}. 

}
\usage{
fundamentalNetworkConcepts(adj, GS = NULL)
}
\arguments{
  \item{adj}{ an adjacency matrix, that is a square, symmetric matrix with entries between 0 and 1}
  \item{GS}{ a node significance measure: a vector of the same length as the number of rows
(and columns) of the adjacency matrix.  }
}
\value{
  A list with the following components:
   \item{Connectivity}{a numerical vector that reports the connectivity (also known as degree) of each
node. This fundamental network concept is also known as whole network connectivity. One can also define
the scaled connectivity \code{K=Connectivity/max(Connectivity)} which is used for computing the hub gene
significance.}
 
  \item{ScaledConnectivity}{the \code{Connectivity} vector scaled by the highest connectivity in the
network, i.e., \code{Connectivity/max(Connectivity)}.}

  \item{ClusterCoef}{a numerical vector that reports the cluster coefficient for each node. This
fundamental network concept measures the cliquishness of each node.}

  \item{MAR}{a numerical vector that reports the maximum adjacency ratio for each node. \code{MAR[i]}
 equals 1 if all non-zero adjacencies between node \code{i} and the remaining network nodes equal 1. This
fundamental network concept is always 1 for nodes of an unweighted network.  This is a useful measure for
weighted networks since it allows one to determine whether a node has high connectivity because of many
weak connections (small MAR) or because of strong (but few) connections (high MAR), see Horvath and Dong
2008. 
}
 
 \item{Density}{the density of the network. }

 \item{Centralization}{the centralization of the network. }
 \item{Heterogeneity}{the heterogeneity of the network. }
}
\references{ 

Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24 

Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression Network Analysis. PLoS Comput Biol
4(8): e1000117

}

\author{ Steve Horvath }
\seealso{ 
\code{\link{conformityBasedNetworkConcepts}} for calculation of conformity based network concepts
for a network adjacency matrix;

\code{\link{networkConcepts}}, for calculation of conformity based and eigennode based network concepts
for a correlation network.
}
\keyword{ misc }
