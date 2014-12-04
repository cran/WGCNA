\name{moduleMergeUsingKME}
\alias{moduleMergeUsingKME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Merge modules and reassign genes using kME.  }
\description{
This function takes an expression data matrix (and other user-defined parameters), calculates the module membership (kME) values, and adjusts the module assignments, merging modules that are not sufficiently distinct and reassigning modules that were originally assigned suboptimally.
}
\usage{
moduleMergeUsingKME(
   datExpr, colorh, ME = NULL, 
   threshPercent = 50, mergePercent = 25, 
   reassignModules = TRUE, 
   convertGrey = TRUE, 
   omitColors = "grey", 
   reassignScale = 1,  
   threshNumber = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{
An expression data matrix, with samples as rows, genes (or probes) as column.
}
  \item{colorh}{
The color vector (module assignments) corresponding to the columns of datExpr.
}
  \item{ME}{
Either NULL (default), at which point the module eigengenes will be calculated, or pre-calculated module eigengenes for each of the modules, with samples as rows (corresponding to datExpr), and modules corresponding to columns (column names MUST be module colors or module colors prefixed by "ME" or "PC").
}
  \item{threshPercent}{
Threshold percent of the number of genes in the module that should be included for the various analyses. For example, in a module with 200 genes, if threshPercent=50 (default), then 50 genes will be checked for reassignment and used to test whether two modules should be merged.  See also threshNumber.
}
  \item{mergePercent}{
If greater than this percent of the assigned genes are above the threshold are in a module other than the assigned module, then these two modules will be merged.  For example, if mergePercent=25 (default), and the 70 out of 200 genes in the blue module were more highly correlated with the black module eigengene, then all genes in the blue module would be reassigned to the black module.
}
  \item{reassignModules}{
If TRUE (default), genes are resassigned to the module with which they have the highest module membership (kME), but only if their kME is above the threshPercent (or threshNumber) threshold of that module.
}
  \item{convertGrey}{
If TRUE (default), unassigned (grey) genes are assigned as in "reassignModules"
}
  \item{omitColors}{
These are all of the module assignments which indicate genes that are not assigned to modules (default="grey").  These genes will all be assigned as "grey" by this function.
}
  \item{reassignScale}{
A value between 0 and 1 (default) which determines how the threshPercent gets scaled for reassigning genes.  Smaller values reassign more genes, but does not affect the merging process.
}
  \item{threshNumber}{
Either NULL (default) or, if entered, every module is counted as having exactly threshNumber genes, and threshPercent it ignored.  This parameter should have the effect of 
}
}
\value{
  \item{moduleColors}{
The NEW color vector (module assignments) corresponding to the columns of datExpr, after module merging and reassignments.
}
  \item{mergeLog}{
A log of the order in which modules were merged, for reference.
}
}
\author{
Jeremy Miller
}
\note{
Note that this function should be considered "experimental" as it has only been beta tested.  Please e-mail jeremyinla@gmail.com if you have any issues with the function.
}
\examples{

## First simulate some data and the resulting network dendrogram
set.seed(100)
MEturquoise = sample(1:100,50)
MEblue      = sample(1:100,50)
MEbrown     = sample(1:100,50)
MEyellow    = sample(1:100,50) 
MEgreen     = c(MEyellow[1:30], sample(1:100,20))
MEred	    = c(MEbrown [1:20], sample(1:100,30))
#MEblack	    = c(MEblue  [1:25], sample(1:100,25))
ME   = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred)#, MEblack)
dat1 = simulateDatExpr(ME, 400, c(0.15,0.13,0.12,0.10,0.09,0.09,0.1), signed=TRUE)
TOM1  = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
tree1 = fastcluster::hclust(as.dist(1-TOM1),method="average")

## Here is an example using different mergePercentages, 
# setting an inclusive threshPercent (91)
colorh1  <- colorPlot <- labels2colors(dat1$allLabels)
merges = c(65,40,20,5)
for (m in merges)  
   colorPlot = cbind(colorPlot, 
                     moduleMergeUsingKME(dat1$datExpr,colorh1,
                        threshPercent=91, mergePercent=m)$moduleColors)
plotDendroAndColors(tree1, colorPlot, c("ORIG",merges), dendroLabels=FALSE)

## Here is an example using a lower reassignScale (so that more genes get reassigned)
colorh1  <- colorPlot <- labels2colors(dat1$allLabels)
merges = c(65,40,20,5)
for (m in merges)  colorPlot = cbind(colorPlot, 
  moduleMergeUsingKME(dat1$datExpr,colorh1,threshPercent=91, 
                      reassignScale=0.7, mergePercent=m)$moduleColors)
plotDendroAndColors(tree1, colorPlot, c("ORIG",merges), dendroLabels=FALSE)

## Here is an example using a less-inclusive threshPercent (75), 
# little if anything is merged.

colorh1  <- colorPlot <- labels2colors(dat1$allLabels)
merges = c(65,40,20,5)
for (m in merges)  colorPlot = cbind(colorPlot, 
  moduleMergeUsingKME(dat1$datExpr,colorh1,
                      threshPercent=75, mergePercent=m)$moduleColors)
plotDendroAndColors(tree1, colorPlot, c("ORIG",merges), dendroLabels=FALSE)
# (Note that with real data, the default threshPercent=50 usually results 
# in some modules being merged)

}
\keyword{misc}
