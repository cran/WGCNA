\name{consensusDissTOMandTree}
\alias{consensusDissTOMandTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Consensus clustering based on topological overlap and hierarchical clustering }
\description{
This function makes a consensus network using all of the default values in the WGCNA library.  Details regarding how consensus modules are formed can be found here: http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.pdf
}
\usage{
consensusDissTOMandTree(multiExpr, softPower, TOM = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{
Expression data in the multi-set format (see checkSets). A vector of lists, one per set. Each set must contain a component data that contains the expression data.  Rows correspond to samples and columns to genes or probes. Two or more sets of data must be included and adjacencies cannot be used.
}
  \item{softPower}{
Soft thresholding power used to make each of the networks in multiExpr.
}
  \item{TOM}{
A LIST of matrices holding the topological overlap corresponding to the sets in multiExpr, if they have already been calculated. Otherwise, keep TOM set as NULL (default), and TOM similarities will be calculated using the WGCNA defaults.  If inputted, this variable must be a list with each entree a TOM corresponding to the same entries in multiExpr.
}
}
\value{
  \item{consensusTOM}{
The TOM difference matrix (1-TOM similarity) corresponding to the consensus network.
}
  \item{consTree}{
Returned value is the same as that of hclust: An object of class hclust which describes the tree produced by the clustering process.  This tree corresponds to the dissimilarity matrix consensusTOM.
}
}
\references{
Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between co-expression modules. BMC Systems Biology 2007, 1:54
}
\author{
Peter Langfelder, Steve Horvath, Jeremy Miller
}
\seealso{
  \code{\link{blockwiseConsensusModules}}
}
\examples{
# Example consensus network using two simulated data sets

set.seed    = 100
MEturquoise = sample(1:100,50)
MEblue      = sample(1:100,50)
MEbrown     = sample(1:100,50)
MEyellow    = sample(1:100,50) 
MEgreen     = sample(1:100,50)

ME   = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen)
dat1 = simulateDatExpr(ME,300,c(0.2,  0.10,  0.10,  0.10,  0.10,  0.2), signed=TRUE)
dat2 = simulateDatExpr(ME,300,c(0.18, 0.11, 0.11, 0.09, 0.11, 0.23),signed=TRUE)
multiExpr = list(S1=list(data=dat1$datExpr),S2=list(data=dat2$datExpr))
softPower=8

consensusNetwork = consensusDissTOMandTree(multiExpr, softPower)
plotDendroAndColors(consensusNetwork$consTree, cbind(labels2colors(dat1$allLabels), 
     labels2colors(dat2$allLabels)),c("S1","S2"), dendroLabels=FALSE)

}
\keyword{misc}
