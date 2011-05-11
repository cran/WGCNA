\name{conformityDecomposition}
\alias{conformityDecomposition}
\title{
Conformity and module based decomposition of a network adjacency matrix.
}
\description{
The function calculates the conformity based approximation \code{A.CF} of an adjacency matrix and a factorizability measure code{Factorizability}. If a module assignment \code{Cl} is provided, it also estimates a corresponding intermodular adjacency matrix. In this case, function automatically carries out the module- and conformity based decomposition of the adjacency matrix described in chapter 2 of (Horvath 2011). }
\usage{
conformityDecomposition(adj, Cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj}{a symmetric numeric matrix (or data frame) whose entries lie between 0 and 1.
}
  \item{Cl}{a vector (or factor variable) of length equal to the number of rows of \code{adj}. The variable assigns each network node (row of \code{adj}) to a module. The entries of \code{Cl} could be integers or character strings.
}
}
\details{We distinguish two situation depending on whether or not \code{Cl} equals \code{NULL}. 
1) Let us start out assuming that \code{Cl = NULL}. In this case,
the function calculates the conformity vector for a general, possibly non-factorizable network \code{adj} by minimizing a quadratic (sums of squares) loss function. The conformity and factorizability for an adjacency matrix is defined in (Dong and Horvath 2007, Horvath and Dong 2008) but we briefly describe it in the following. A network is called exactly factorizable if the pairwise connection strength (adjacency) between 2 network nodes can be factored into node specific contributions, named node 'conformity', i.e. if \code{adj[i,j]=Conformity[i]*Conformity[j]}. The conformity turns out to be highly related to the network connectivity (aka degree). If \code{adj} is not exactly factorizable, then the function \code{conformityDecomposition} calculates a conformity vector of the exactly factorizable network that best approximates \code{adj}. The factorizability measure \code{Factorizability} is a number between 0 and 1. The higher \code{Factorizability}, the more factorizable is \code{adj}. Warning: the algorithm may only converge to a local optimum and it may not converge at all. Also see the notes below.

2) Let us now assume that \code{Cl} is not NULL, i.e. it specifies the module assignment of each node.
Then the function calculates a module- and CF-based approximation of \code{adj} (explained in chapter 2 in Horvath 2011). In this case, the function calculates a conformity vector \code{Conformity} and a matrix \code{IntermodularAdjacency} such that \code{adj[i,j]} is approximately equal to
\code{Conformity[i]*Conformity[j]*IntermodularAdjacency[module.index[i],module.index[j]]} where \code{module.index[i]} is the row of the matrix \code{IntermodularAdjacency} that corresponds to the module assigned to node i.
To estimate \code{Conformity} and a matrix \code{IntermodularAdjacency}, the function attempts to minimize a quadratic loss function (sums of squares).
Currently, the function only implements a heuristic algorithm for optimizing the objective function (chapter 2 of Horvath 2011). Another, more accurate Majorization Minorization (MM) algorithm for the decomposition is implemented in the function \code{propensityDecomposition} by Ranola et al (2011).
}
\value{
\item{A.CF}{a symmetric matrix that approximates the input matrix \code{adj}. Roughly speaking, the i,j-the element of the matrix equals \code{Conformity[i]*Conformity[j]*IntermodularAdjacency[module.index[i],module.index[j]]} where \code{module.index[i]} is the row of the matrix \code{IntermodularAdjacency} that corresponds to the module assigned to node i. } 
\item{Conformity}{a numeric vector whose entries correspond to the rows of code{adj}. If \code{Cl=NULL} then \code{Conformity[i]} is the conformity.
If \code{Cl} is not NULL then \code{Conformity[i]} is the intramodular conformity with respect to the module that node i belongs to.
}
\item{IntermodularAdjacency}{ a symmetric matrix (data frame) whose rows and columns correspond to the number of modules specified in \code{Cl}. Interpretation: it measures the similarity (adjacency) between the modules. In this case, the rows (and columns) of \code{IntermodularAdjacency} correspond to the entries of \code{Cl.level}. }
\item{Factorizability}{ is a number between 0 and 1. If \code{Cl=NULL} then it equals 1, if (and only if) \code{adj} is exactly factorizable. If \code{Cl} is a vector, then it measures how well the module- and CF based decomposition approximates \code{adj}.  }
\item{Cl.level}{ is a vector of character strings which correspond to the factor levels of the module assignment \code{Cl}. Incidentally, the function automatically turns \code{Cl} into a factor variable. The components of Conformity and \code{IntramodularFactorizability} correspond to the entries of \code{Cl.level}. }
\item{IntramodularFactorizability}{ is a numeric vector of length equal to the number of modules specified by \code{Cl}. Its entries report the factorizability measure for each module. The components correspond to the entries of \code{Cl.level}.}
\item{listConformity}{} 
}
\references{
Dong J, Horvath S (2007) Understanding Network Concepts in Modules. BMC Systems Biology 2007, June 1:24
Horvath S, Dong J (2008) Geometric Interpretation of Gene Co-Expression Network Analysis. PloS Computational Biology. 4(8): e1000117. PMID: 18704157
Horvath S (2011) Weighted Network Analysis. Applications in Genomics and Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8
Ranola JMO, Langfelder P, Song L, Horvath S, Lange K (2011) An MM algorithm for the module- and propensity based decomposition of a network. Currently a draft.
 }
\author{
Steve Horvath
}
\note{Regarding the situation when \code{Cl=NULL}.
One can easily show that the conformity vector is not unique if \code{adj} contains only 2 nodes. However, for more than 2 nodes the conformity is uniquely defined when dealing with an exactly factorizable weighted network whose entries \code{adj[i,j]} are larger than 0. In this case, one can get explicit formulas for the conformity (Dong and Horvath 2007). 
}

\seealso{
\code{\link{conformityBasedNetworkConcepts}}
%%\code{\link{ propensityDecomposition }}
}
\examples{

# assume the number of nodes can be divided by 2 and by 3
n=6
# here is a perfectly factorizable matrix
A=matrix(1,nrow=n,ncol=n)
# this provides the conformity vector and factorizability measure
conformityDecomposition(adj=A)
# now assume we have a class assignment
Cl=rep(c(1,2),c(n/2,n/2))
conformityDecomposition(adj=A,Cl=Cl)
# here is a block diagonal matrix
blockdiag.A=A
blockdiag.A[1:(n/3),(n/3+1):n]=0
blockdiag.A[(n/3+1):n , 1:(n/3)]=0
block.Cl=rep(c(1,2),c(n/3,2*n/3))
conformityDecomposition(adj= blockdiag.A,Cl=block.Cl)

# another block diagonal matrix
blockdiag.A=A
blockdiag.A[1:(n/3),(n/3+1):n]=0.3
blockdiag.A[(n/3+1):n , 1:(n/3)]=0.3
block.Cl=rep(c(1,2),c(n/3,2*n/3))
conformityDecomposition(adj= blockdiag.A,Cl=block.Cl)

}
\keyword{misc }
