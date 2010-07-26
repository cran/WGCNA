\name{networkConcepts}
\alias{networkConcepts}
\title{ Calculations of network concepts}
\description{
  This functions calculates various network concepts (topological properties, network indices) of a
network calculated from expression data. See details for a detailed description.
}
\usage{
networkConcepts(datExpr, power = 1, trait = NULL, networkType = "unsigned")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{ a data frame containg the expression data, with rows corresponding to samples and
columns to genes (nodes). }

  \item{power}{ soft thresholding power.}
  \item{trait}{optional specification of a sample trait. A vector of length equal the number of samples
in \code{datExpr}. }
  \item{networkType}{ network type. Recognized values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, and \code{"signed hybrid"}. }
}


\details{
This function computes various network concepts (also known as network statistics, topological
properties, or network indices) for a weighted correlation network. The nodes of the weighted correlation
network will be constructed between the columns (interpreted as nodes) of the input \code{datExpr}. 
If the option
\code{networkType="unsigned"} then the adjacency between nodes i and j is defined as
\code{A[i,j]=abs(cor(datExpr[,i],datExpr[,j]))^power}.
In the following, we use the term gene and node interchangeably since these methods were originally
developed for gene networks. The function computes the following
4 types of network concepts (introduced in Horvath and Dong 2008):

Type I: fundamental network concepts are defined as a function of the off-diagonal elements of an
adjacency matrix A and/or a node significance measure GS. These network concepts can be defined for  any
network (not just correlation networks).
The adjacency matrix of an unsigned weighted correlation network is given by
\code{A=abs(cor(datExpr,use="p"))^power} and the trait based gene significance measure is given by 
\code{GS= abs(cor(datExpr,trait, use="p"))^power} where \code{datExpr}, \code{trait}, \code{power}
 are input parameters.

Type II: conformity-based network concepts are functions of the off-diagonal elements of the conformity
based adjacency matrix \code{A.CF=CF*t(CF)} and/or the node significance measure. These network concepts
are
defined for any network for which a conformity vector can be defined. Details: For any adjacency matrix
\code{A}, the conformity vector \code{CF} is calculated by requiring that \code{A[i,j]} is 
approximately equal to \code{CF[i]*CF[j]}.
Using the conformity one can define the matrix \code{A.CF=CF*t(CF)} which is the outer product of 
the conformity
vector with itself. In general, \code{A.CF} is not an adjacency matrix since its diagonal elements 
are different
from 1. If the off-diagonal elements of \code{A.CF} are similar to those of \code{A}
according to the Frobenius matrix
norm, then \code{A} is approximately factorizable. To measure the factorizability of a network, one can
calculate the \code{Factorizability}, which is a number between 0 and 1 (Dong and Horvath 2007). T
he conformity
is defined using a monotonic, iterative algorithm that maximizes the factorizability measure. 

Type III: approximate conformity based network concepts are functions of all elements of the conformity
based adjacency matrix \code{A.CF} (including the diagonal) and/or the node significance measure
\code{GS}. These
network concepts are very useful for deriving relationships between network concepts in networks that are
approximately factorizable.

Type IV: eigengene-based (also known as eigennode-based) network concepts are functions of the
eigengene-based adjacency matrix \code{A.E=ConformityE*t(ConformityE)} (diagonal included) and/or the
corresponding eigengene-based gene significance measure \code{GSE}.  These network concepts can only be 
defined
for correlation networks. Details: The columns (nodes) of \code{datExpr} can be summarized with the 
first principal
component, which is referred to as Eigengene in coexpression network analysis. In general correlation
networks, it is called eigennode. The eigengene-based conformity \code{ConformityE[i]} is defined as
\code{abs(cor(datE[,i], Eigengene))^power} where the power corresponds to the power used for defining the
weighted adjacency matrix \code{A}. The eigengene-based conformity can also be used to define an 
eigengene-based
adjacency matrix \code{A.E=ConformityE*t(ConformityE)}. 
The eigengene based factorizability \code{EF(datE)} is a number between 0 and 1 that measures how well
\code{A.E}
approximates \code{A} when the power parameter equals 1. \code{EF(datE)} is defined with respect to the
singular values
of \code{datExpr}. For a trait based node significance measure \code{GS=abs(cor(datE,trait))^power}, 
one can also define
an eigengene-based node significance measure \code{GSE[i]=ConformityE[i]*EigengeneSignificance} where the
eigengene significance \code{abs(cor(Eigengene,trait))^power} is defined as power of the absolute value
of the
correlation between eigengene and trait.
Eigengene-based network concepts are very useful for providing a geometric interpretation of network
concepts and for deriving relationships between network concepts. For example, the hub gene significance
measure and its eigengene-based analog have been used to characterize networks where highly connected hub
genes are important with regard to a trait based gene significance measure (Horvath and Dong 2008).

}
\value{
  A list with the following components:
  \item{Summary}{a data frame whose rows report network concepts that only depend on the adjacency
matrix.
Density (mean adjacency), Centralization , Heterogeneity (coefficient of variation of the connectivity),
Mean ClusterCoef, Mean Connectivity.
The columns of the data frame report the 4 types of network concepts mentioned in the description:
Fundamental concepts, eigengene-based concepts, conformity-based concepts, and approximate
conformity-based concepts.}

  \item{Size}{reports the network size, i.e. the number of nodes, which equals the number of columns of
the input data frame \code{datExpr}.}

  \item{Factorizability}{a number between 0 and 1. The closer it is to 1, the better the off-diagonal
elements of the conformity based network \code{A.CF} approximate those of \code{A} 
(according to the Frobenius norm). }

  \item{Eigengene}{the first principal component of the standardized columns of \code{datExpr}. The
number of
components of this vector equals the number of rows of \code{datExpr}.}

  \item{VarExplained}{the proportion of variance explained by the first principal component (the
\code{Eigengene}). It is numerically different from the eigengene based factorizability. 
While \code{VarExplained} is
based on the squares of the singular values of \code{datExpr}, 
the eigengene-based factorizability is based on
fourth powers of the singular values. }

  \item{Conformity}{numerical vector giving the conformity. 
The number of components of the conformity vector equals the number of columns in
\code{datExpr}. The conformity is often highly correlated with the vector of node connectivities. The
conformity is computed using an iterative algorithm for maximizing the factorizability measure. The
algorithm and related network concepts are described in Dong and Horvath 2007.} 

  \item{ClusterCoef}{a numerical vector that reports the cluster coefficient for each node. This
fundamental network concept measures the cliquishness of each node.}

  \item{Connectivity}{a numerical vector that reports the connectivity (also known as degree) of each
node. This fundamental network concept is also known as whole network connectivity. One can also define
the scaled connectivity \code{K=Connectivity/max(Connectivity)} which is used for computing the hub gene
significance.}

  \item{MAR}{a numerical vector that reports the maximum adjacency ratio for each node. \code{MAR[i]}
 equals 1
if all non-zero adjacencies between node \code{i} and the remaining network nodes equal 1. This
fundamental
network concept is always 1 for nodes of an unweighted network.  This is a useful measure for weighted
networks since it allows one to determine whether a node has high connectivity because of many weak
connections (small MAR) or because of strong (but few) connections (high MAR), see Horvath and Dong 2008.
}

  \item{ConformityE}{a numerical vector that reports the eigengene based (aka eigenenode based)
conformity for the correlation network. The number of components equals the number of columns of
\code{datExpr}.}

  \item{GS}{a numerical vector that encodes the node (gene) significance. The i-th component equals the
node significance of the i-th column of \code{datExpr} if a sample trait was supplied to the function
(input
trait). \code{GS[i]=abs(cor(datE[,i], trait, use="p"))^power} }

  \item{GSE}{a numerical vector that reports the eigengene based gene significance measure.  Its i-th
component is given by \code{GSE[i]=ConformityE[i]*EigengeneSignificance} where the eigengene significance
\code{abs(cor(Eigengene,trait))^power} is defined as power of the absolute value of the correlation
between
eigengene and trait.}

  \item{Significance}{a data frame whose rows report network concepts that also depend on the trait based
node significance measure. The rows correspond to network concepts and the columns correspond to the type
of network concept (fundamental versus eigengene based). The first row of the data frame reports the
network significance. The fundamental version of this network concepts is the average gene
significance=mean(GS). The eigengene based analog of this concept is defined as mean(GSE). The second row
reports the hub gene significance which is defined as slope of the intercept only regression model that
regresses the gene significance on the scaled network connectivity K. The third row reports the eigengene
significance \code{abs(cor(Eigengene,trait))^power}. More details can be found in Horvath and Dong
(2008).}

}
\references{
Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17

Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24

Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression Network Analysis. PLoS Comput Biol
4(8): e1000117
}

\seealso{ 

  \code{\link{conformityBasedNetworkConcepts}} for approximate conformity-based network concepts

  \code{\link{fundamentalNetworkConcepts}} for calculation of fundamental network concepts only.

}

\author{ Jun Dong, Steve Horvath, Peter Langfelder }

\keyword{ misc }

