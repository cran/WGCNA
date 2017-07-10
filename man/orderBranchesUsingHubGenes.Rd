\name{orderBranchesUsingHubGenes}
\alias{orderBranchesUsingHubGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Optimize dendrogram using branch swaps and reflections.  }
\description{
This function takes as input the hierarchical clustering tree as well as a subset of genes in the network (generally corresponding to branches in the tree), then returns a semi-optimally ordered tree.  The idea is to maximize the correlations between adjacent branches in the dendrogram, in as much as that is possible by adjusting the arbitrary positionings of the branches by swapping and reflecting branches.
}
\usage{
orderBranchesUsingHubGenes(
  hierTOM, 
  datExpr = NULL, colorh = NULL, 
  type = "signed", adj = NULL, iter = NULL, 
  useReflections = FALSE, allowNonoptimalSwaps = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hierTOM}{
A hierarchical clustering object (or gene tree) that is used to plot the dendrogram.  For example, the output object from the function hclust or fastcluster::hclust.  Note that elements of hierTOM$order MUST be named (for example, with the corresponding gene name).
}
  \item{datExpr}{
Gene expression data with rows as samples and columns as genes, or NULL if a pre-made adjacency is entered.  Column names of datExpr must be a subset of gene names of hierTOM$order.
}
  \item{colorh}{
The module assignments (color vectors) corresponding to the rows in datExpr, or NULL if a pre-made adjacency is entered.
}
  \item{type}{
What type of network is being entered.  Common choices are "signed" (default) and "unsigned".  With "signed" negative correlations count against, whereas with "unsigned" negative correlations are treated identically as positive correlations.
}
  \item{adj}{
Either NULL (default) or an adjacency (or any other square) matrix with rows and columns corresponding to a subset of the genes in hierTOM$order.  If entered, datExpr, colorh, and type are all ignored.  Typically, this would be left blank but could include correlations between module eigengenes, with rows and columns renamed as genes in the corresponding modules, for example.
}
  \item{iter}{
The number of iterations to run the function in search of optimal branch ordering.  The default is the square of the number of modules (or the quare of the number of genes in the adjacency matrix).
}
  \item{useReflections}{
If TRUE, both reflections and branch swapping will be used to optimize dendrogram.  If FALSE (default) only branch swapping will be used. 
}
  \item{allowNonoptimalSwaps}{
If TRUE, there is chance (that decreases with each iteration) of swapping / reflecting branches whether or not the new correlation between expression of genes in adjacent branches is better or worse.  The idea (which has not been sufficiently tested), is that this would prevent the function from getting stuck at a local maxima of correlation.  If FALSE (default), the swapping / reflection of branches only occurs if it results in a higher correlation between adjacent branches.
}
}
\value{
  \item{hierTOM}{
A hierarchical clustering object with the hierTOM$order variable properly adjusted, but all other variables identical as the heirTOM input.
}
  \item{changeLog}{
A log of all of the changes that were made to the dendrogram, including what change was made, on what iteration, and the Old and New scores based on correlation.  These scores have arbitrary units, but higher is better.
}
}
\author{
Jeremy Miller
}
\note{
This function is very slow and is still in an *experimental* function.  We have not had problems with ~10 modules across ~5000 genes, although theoretically it should work for many more genes and modules, depending upon the speed of the computer running R.  Please address any problems or suggestions to jeremyinla@gmail.com.
}
\examples{

\dontrun{
## Example: first simulate some data.

MEturquoise = sample(1:100,50)
MEblue      = c(MEturquoise[1:25], sample(1:100,25))
MEbrown     = sample(1:100,50)
MEyellow    = sample(1:100,50) 
MEgreen     = c(MEyellow[1:30], sample(1:100,20))
MEred	    = c(MEbrown [1:20], sample(1:100,30))
ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred)
dat1   = simulateDatExpr(ME,400,c(0.16,0.12,0.11,0.10,0.10,0.10,0.1), signed=TRUE)
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
# and output pdf with resulting images

pdf("DENDROGRAM_PLOTS.pdf",width=10,height=5)
plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Starting Dendrogram")

tree1 = swapTwoBranches(tree1,hubs["red"],hubs["turquoise"])
plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Swap blue/turquoise and red/brown")

tree1 = reflectBranch(tree1,hubs["blue"],hubs["green"])
plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Reflect turquoise/blue")

# (This function will take a few minutes)
out = orderBranchesUsingHubGenes(tree1,datExpr,colorh2,useReflections=TRUE,iter=100)
tree1 = out$geneTree
plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Semi-optimal branch order")

out$changeLog

dev.off()
}
}

\keyword{misc}
