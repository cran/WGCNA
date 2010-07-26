\name{conformityBasedNetworkConcepts}
\alias{conformityBasedNetworkConcepts}
\title{ Calculation of conformity-based network concepts. }
\description{
This function computes 3 types of network concepts (also known as network indices or statistics) based on
an adjacency matrix and optionally a node significance measure.
}
\usage{
conformityBasedNetworkConcepts(adj, GS = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj}{ adjacency matrix. A symmetric matrix with components between 0 and 1. }
  \item{GS}{ optional node significance measure. A vector with length equal the dimension of \code{adj}. }
}
\details{
This function computes 3 types of network concepts (also known as network indices or statistics) based on
an adjacency matrix and optionally a node significance measure. Specifically, it computes I) fundamental
network concepts, II) conformity based network concepts, and III) approximate conformity based network
concepts.
These network concepts are defined for any symmetric adjacency matrix (weighted and unweighted). The
network concepts are described in Dong and Horvath (2007) and Horvath and Dong (2008).
In the following, we use the term gene and node interchangeably since these methods were originally
developed for gene networks. In the following, we briefly describe the 
3 types of network concepts:

Type I: fundamental network concepts are defined as a function of the off-diagonal elements of an
adjacency matrix A and/or a node significance measure GS. 
Type II: conformity-based network concepts are functions of the off-diagonal elements of the conformity
based adjacency matrix A.CF=CF*t(CF) and/or the node significance measure. These network concepts are
defined for any network for which a conformity vector can be defined. Details: For any adjacency matrix
A, the conformity vector CF is calculated by requiring that A[i,j] is approximately equal to CF[i]*CF[j].
Using the conformity one can define the matrix A.CF=CF*t(CF) which is the outer product of the conformity
vector with itself. In general, A.CF is not an adjacency matrix since its diagonal elements are different
from 1. If the off-diagonal elements of A.CF are similar to those of A according to the Frobenius matrix
norm, then A is approximately factorizable. To measure the factorizability of a network, one can
calculate the Factorizability, which is a number between 0 and 1 (Dong and Horvath 2007). The conformity
is defined using a monotonic, iterative algorithm that maximizes the factorizability measure. 
Type III: approximate conformity based network concepts are functions of all elements of the conformity
based adjacency matrix A.CF (including the diagonal) and/or the node significance measure GS. These
network concepts are very useful for deriving relationships between network concepts in networks that are
approximately factorizable.

}
\value{
  A list with the following components:

  \item{Factorizability}{number between 0 and 1 giving the factorizability of the matrix. 
The closer to 1 the higher the evidence of factorizability,
that is, A-I is close to outer(CF,CF)-diag(CF^2).}

  \item{fundamentalNCs}{fundamental network concepts, that is network concepts calculated directly from
the given adjacency matrix \code{adj}. A list with components \code{ScaledConnectivity} (giving the
scaled connectivity of each node), \code{Connectivity} (connectivity of each node), \code{ClusterCoef}
(the clustering coefficient of each node), \code{MAR} (maximum adjacency ratio of each node),
\code{Density} (the mean density of the network), \code{Centralization} (the centralization of the
network), \code{Heterogeneity} (the heterogeneity of the network). If the input node significance
\code{GS} is specified, the following additional components are included: \code{NetworkSignificance}
(network significance, the mean node significance), and \code{HubNodeSignificance} (hub node significance
given by the linear regression of node significance on connectivity). }

  \item{conformityBasedNCs}{network concepts based on an approximate adjacency matrix given by the
outer product of the conformity vector but with unit diagonal. A list with components \code{Conformity}
(the conformity vector) and \code{Connectivity.CF, ClusterCoef.CF, MAR.CF, Density.CF, Centralization.CF,
Heterogeneity.CF} giving the conformity-based analogs of the above network concepts. }

  \item{approximateConformityBasedNCs}{network concepts based on an approximate adjacency matrix given by
the outer product of the conformity vector. A list with components \code{Conformity} 
(the conformity vector) and \code{Connectivity.CF.App, ClusterCoef.CF.App, MAR.CF.App, Density.CF.App,
Centralization.CF.App,
Heterogeneity.CF.App} giving the conformity-based analogs of the above network concepts. }

}
\references{ 
  Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24 
Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression Network Analysis. PLoS Comput Biol
4(8): e1000117 
}

\author{ Steve Horvath }
\seealso{
   \code{\link{networkConcepts}} for calculation of eigennode based network concepts for a
   correlation network;

   \code{\link{fundamentalNetworkConcepts}} for calculation of fundamental network concepts only. 
}

\keyword{ misc }% __ONLY ONE__ keyword per line
