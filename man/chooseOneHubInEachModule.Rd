\name{chooseOneHubInEachModule}
\alias{chooseOneHubInEachModule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Chooses a single hub gene in each module }
\description{
chooseOneHubInEachModule returns one gene in each module with high connectivity, given a number of randomly selected genes to test.
}
\usage{
chooseOneHubInEachModule(
   datExpr, 
   colorh,  
   numGenes = 100, 
   omitColors = "grey", 
   power = 2, 
   type = "signed", 
   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{
Gene expression data with rows as samples and columns as genes.
}
  \item{colorh}{
The module assignments (color vectors) corresponding to the rows in datExpr.
}
  \item{numGenes}{
Th number of random genes to select per module.  Higher number of genes increases the accuracy of hub selection but slows down the function.
}
  \item{omitColors}{
All colors in this character vector (default is "grey") are ignored by this function.
}
  \item{power}{
Power to use for the adjacency network (default = 2).
}
  \item{type}{
What type of network is being entered.  Common choices are "signed" (default) and "unsigned".  With "signed" negative correlations count against, whereas with "unsigned" negative correlations are treated identically as positive correlations.
}
  \item{\dots}{
Any other parameters accepted by the *adjacency* function
}
}
\value{
Both functions output a character vector of genes, where the genes are the hub gene picked for each module, and the names correspond to the module in which each gene is a hub.
}
\author{
Jeremy Miller
}
\examples{
## Example: first simulate some data.

MEturquoise = sample(1:100,50)
MEblue      = sample(1:100,50)
MEbrown     = sample(1:100,50)
MEyellow    = sample(1:100,50) 
MEgreen     = c(MEyellow[1:30], sample(1:100,20))
MEred	    = c(MEbrown [1:20], sample(1:100,30))
MEblack	    = c(MEblue  [1:25], sample(1:100,25))
ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred, MEblack)
dat1   = simulateDatExpr(ME,300,c(0.2,0.1,0.08,0.051,0.05,0.042,0.041,0.3), 
                         signed=TRUE)
TOM1   = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
colnames(TOM1) <- rownames(TOM1) <- colnames(dat1$datExpr)
tree1 <- tree2 <- fastcluster::hclust(as.dist(1-TOM1),method="average")
colorh = labels2colors(dat1$allLabels)
hubs    = chooseOneHubInEachModule(dat1$datExpr, colorh)
hubs

}

\keyword{misc}
