\name{swapTwoBranches}
\alias{swapTwoBranches}
\alias{reflectBranch}
\alias{selectBranch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Select, swap, or reflect branches in a dendrogram.  }
\description{
swapTwoBranches takes the a gene tree object and two genes as input, and swaps the branches containing these two genes at the nearest branch point of the dendrogram.

reflectBranch takes the a gene tree object and two genes as input, and reflects the branch containing the first gene at the nearest branch point of the dendrogram.

selectBranch takes the a gene tree object and two genes as input, and outputs indices for all genes in the branch containing the first gene, up to the nearest branch point of the dendrogram.
}
\usage{
swapTwoBranches(hierTOM, g1, g2)
reflectBranch(hierTOM, g1, g2, both = FALSE)
selectBranch(hierTOM, g1, g2)
}
\arguments{
  \item{hierTOM}{
A hierarchical clustering object (or gene tree) that is used to plot the dendrogram.  For example, the output object from the function hclust or fastcluster::hclust.  Note that elements of hierTOM$order MUST be named (for example, with the corresponding gene name).
}
  \item{g1}{
Any gene in the branch of interest.
}
  \item{g2}{
Any gene in a branch directly adjacent to the branch of interest.
}
  \item{both}{
Logical: should the selection include the branch gene \code{g2}?
}
}

\value{
swapTwoBranches and reflectBranch return a hierarchical clustering object with the hierTOM$order variable properly adjusted, but all other variables identical as the heirTOM input.

selectBranch returns a numeric vector corresponding to all genes in the requested branch.
}
\author{
Jeremy Miller
}
\examples{
\dontrun{
## Example: first simulate some data.
n = 30;
n2 = 2*n;
n.3 = 20;
n.5 = 10;
MEturquoise = sample(1:(2*n),n)
MEblue      = c(MEturquoise[1:(n/2)], sample(1:(2*n),n/2))
MEbrown     = sample(1:n2,n)
MEyellow    = sample(1:n2,n) 
MEgreen     = c(MEyellow[1:n.3], sample(1:n2,n.5))
MEred	    = c(MEbrown [1:n.5], sample(1:n2,n.3))

ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred)
dat1   = simulateDatExpr(ME,8*n ,c(0.16,0.12,0.11,0.10,0.10,0.09,0.15), 
                         signed=TRUE)
TOM1   = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
colnames(TOM1) <- rownames(TOM1) <- colnames(dat1$datExpr)
tree1  = fastcluster::hclust(as.dist(1-TOM1),method="average")
colorh = labels2colors(dat1$allLabels)

plotDendroAndColors(tree1,colorh,dendroLabels=FALSE)

## Reassign modules using the selectBranch and chooseOneHubInEachModule functions

datExpr = dat1$datExpr
hubs    = chooseOneHubInEachModule(datExpr, colorh)
colorh2 = rep("grey", length(colorh))
colorh2 [selectBranch(tree1,hubs["blue"],hubs["turquoise"])] = "blue"
colorh2 [selectBranch(tree1,hubs["turquoise"],hubs["blue"])] = "turquoise"
colorh2 [selectBranch(tree1,hubs["green"],hubs["yellow"])]   = "green"
colorh2 [selectBranch(tree1,hubs["yellow"],hubs["green"])]   = "yellow"
colorh2 [selectBranch(tree1,hubs["red"],hubs["brown"])]      = "red"
colorh2 [selectBranch(tree1,hubs["brown"],hubs["red"])]      = "brown"
plotDendroAndColors(tree1,cbind(colorh,colorh2),c("Old","New"),dendroLabels=FALSE)

## Now swap and reflect some branches, then optimize the order of the branches

# Open a suitably sized graphics window

sizeGrWindow(12,9);

# partition the screen for 3 dendrogram + module color plots

layout(matrix(c(1:6), 6, 1), heights = c(0.8, 0.2, 0.8, 0.2, 0.8, 0.2));

plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Starting Dendrogram", 
                    setLayout = FALSE)

tree1 = swapTwoBranches(tree1,hubs["red"],hubs["turquoise"])
plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Swap blue/turquoise and red/brown", 
                    setLayout = FALSE)

tree1 = reflectBranch(tree1,hubs["blue"],hubs["green"])
plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Reflect turquoise/blue", 
                    setLayout = FALSE)

}}
\keyword{misc}
