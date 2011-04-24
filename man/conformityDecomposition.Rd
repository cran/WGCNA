\name{conformityDecomposition}
\alias{conformityDecomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conformity vector(s) and factorizability measure(s) of a network
}
\description{
The function calculates the conformity vector CF and factorizability measure \code{Factorizability}  (number between 0 and 1) for the input adjacency matrix \code{adj}. If a module (class) assignment vector \code{Cl} is provided then the function automatically carries out the module- and conformity based decomposition of the adjacency matrix. In this case, it outputs a list of (intramodular) conformity vectors, an estimate of the intermodular adjacency matrix,
the overall \code{Factorizability} (number between 0 and 1), and a vector of intramodular factorizabilities. 
}
\usage{
conformityDecomposition(adj, Cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj}{ a symmetric numeric matrix (or data frame) whose entries lie between 0 and 1.}
  \item{Cl}{a vector (or factor variable) of length equal to the number of rows of \code{adj}. The variable assigns each network node (row of \code{adj}) to a  a module. The entries of \code{Cl} could be integers or character strings.
}
}
\details{
This function calculates the conformity vector for a general, possibly non-factorizable network by minimizing a quadratic loss function. Warning: the algorithm may only converge to a local optimum and it may not converge at all.
The conformity and factorizability for an adjacency matrix is defined in (Dong and Horvath 2007, Horvath and Dong 2008) but we briefly describe it in the following. A network is called exactly factorizable if the pairwise connection strength (adjacency) between 2 network nodes can be factored into node specific contributions, named node 'conformity', i.e. if \code{adj[i,j]}=CF[i]*CF[j] where CF is the conformity vector. The conformity CF turns out to be highly related to the network connectivity (aka degree). If \code{adj} is not exactly factorizable, then the function \code{conformityDecomposition} calculates a conformity vector of the exactly factorizable network that best approximates \code{adj}. The factorizability measure \code{Factorizability} is a number between 0 and 1. The higher \code{Factorizability}, the more factorizable is \code{adj}. One can easily show that the conformity vector is not unique if \code{adj} contains only 2 nodes. However, for more than 2 nodes the conformity is uniquely defined when dealing with an exactly factorizable weighted network whose entries \code{adj[i,j]} are larger than 0. In this case, one can get explicit formulas for the conformity (Dong and Horvath 2007). 
If the module assignment vector \code{Cl} is provided then the function calculates the module- and CF-based approximation of \code{adj} (explained in chapter 2 in Horvath 2011). Currently, the function only implements a heuristic algorithm for arriving at approximate solutions of the module- and CF-based approximation of \code{adj} (chapter 2 of Horvath 2011). Another, more accurate Majorization Minorization (MM) algorithm for the decomposition is implemented in the function \code{propensityDecomposition} by Ranola et al (2011).

}
\value{
A list with the following components
\item{ Conformity }{If \code{Cl=NULL} then it is a numeric vector of length equal to the number of nodes. Otherwise code{Conformity} is a list whose components provide one (intramodular) conformity vector for each module (class level) specified in \code{Cl}. In this case, the list components 
correspond to the entries of \code{Cl.level}.
 }
\item{ IntermodularAdjacency }{ A symmetric data frame whose rows and columns correspond to the number of modules specified in \code{Cl}. Interpretation: it measures the similarity (adjacency) between the modules. 
In this case, the rows (and columns) of \code{IntermodularAdjacency}
correspond to the entries of \code{Cl.level}.
}
\item{Factorizability }{is a number between 0 and 1. If \code{Cl=NULL} then it equals 1, if (and only if) \code{adj} is exactly factorizable.
If \code{Cl} is a vector, then it measures how well the module- and CF based decomposition approximates \code{adj}.}

\item{IntramodularFactorizability}{is a numeric vector of length equal to the number of modules specified by \code{Cl}. Its entries report the factorizability measure for each module. The components correspond to the entries of \code{Cl.level}.
 }

\item{Cl.level}{ is a vector of character strings which correspond to the factor levels of the module assignment \code{Cl}. Incidentally, the function automatically turns \code{Cl} into a factor variable. The components of \code{Conformity} and
\code{IntramodularFactorizability} correspond to the entries of \code{Cl.level}.
 }
}
\references{
Dong J, Horvath S (2007) Understanding Network Concepts in Modules. BMC Systems Biology 2007, June 1:24

Horvath S, Dong J (2008) Geometric Interpretation of Gene Co-Expression Network Analysis. PloS Computational Biology. 4(8): e1000117. PMID: 18704157

Horvath S (2011) Weighted Network Analysis. Applications in Genomics and Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8

Ranola JMO, Langfelder P, Song L, Horvath S, Lange K (2011) An MM algorithm for the module- and propensity based decomposition of a network. Currently a draft.

%% ~put references to the literature/web site here ~
}
\author{
Steve Horvath
}
\seealso{
\code{conformityBasedNetworkConcepts},\code{propensityDecomposition} 
%\code{propensityClustering}

}
\examples{
# the examples assume that number of nodes is divisible by 6
n=12
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}% __ONLY ONE__ keyword per line
