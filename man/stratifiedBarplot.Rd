\name{stratifiedBarplot}
\alias{stratifiedBarplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bar plots of data across two splitting parameters }
\description{
This function takes an expression matrix which can be split using two separate splitting parameters (ie, control vs AD with multiple brain regions), and plots the results as a barplot. Group average, standard deviations, and relevant Kruskal-Wallis p-values are returned.
}
\usage{
stratifiedBarplot(
  expAll, 
  groups, split, subset, 
  genes = NA, 
  scale = "N", graph = TRUE, 
  las1 = 2, cex1 = 1.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expAll}{
An expression matrix, with rows as samples and genes/probes as columns.  If genes=NA, then column names must be included.
}
  \item{groups}{
A character vector corresponding to the samples in expAll, with each element the group name of the relevant sample or NA for samples not in any group.  For, example: NA, NA, NA, Con, Con, Con, Con, AD, AD, AD, AD, NA, NA.  This trait will be plotted as adjacent bars for each split.
}
  \item{split}{
A character vector corresponding to the samples in expAll, with each element the group splitting name of the relevant sample or NA for samples not in any group.  For, example: NA, NA, NA, Hip, Hip, EC, EC, Hip, Hip, EC, EC, NA, NA.  This trait will be plotted as the same color across each split of the barplot.  For the function to work properly, the same split values should be inputted for each group.
}
  \item{subset}{
A list of one or more genes to compare the expression with.  If the list contains more than one gene, the first element contains the group name. For example, Ribosomes, RPL3, RPL4, RPS3.
}
  \item{genes}{
If entered, this parameter is a list of gene/probe identifiers corresponding to the columns in expAll.
}
  \item{scale}{
For subsets of genes that include more than one gene, this parameter determines how the genes are combined into a single value.  Currently, there are five options: 1) ("N")o scaling (default); 2) first divide each gene by the ("A")verage across samples; 3) first scale genes to ("Z")-score across samples; 4) only take the top ("H")ub gene (ignore all but the highest-connected gene); and 5) take the ("M")odule eigengene.  Note that these scaling methods have not been sufficiently tested, and should be considered experimental.
}
  \item{graph}{
If TRUE (default), bar plot is made.  If FALSE, only the results are returned, and no plot is made.
}
  \item{cex1}{
Sets the graphing parameters of cex.axis and cex.names (default=1.5)
}
  \item{las1}{
Sets the graphing parameter las (default=2).
}
  \item{\dots}{
Other graphing parameters allowed in the barplot function.  Note that the parameters for cex.axis, cex.names, and las are superseded by cex1 and las1 and will therefore be ignored.
}
}
\value{
  \item{splitGroupMeans}{
The group/split averaged expression across each group and split combination.  This is the height of the bars in the graph.
}
  \item{splitGroupSDs}{
The standard deviation of group/split expression across each group and split combination.  This is the height of the error bars in the graph.
}
  \item{splitPvals}{
Kruskal-Wallis p-values for each splitting parameter across groups.
}
  \item{groupPvals}{
Kruskal-Wallis p-values for each group parameter across splits.
}
}
\author{
Jeremy Miller
}
\seealso{
\code{\link{barplot}}, \code{\link{verboseBarplot}}
}
\examples{
# Example: first simulate some data
set.seed(100)
ME.A = sample(1:100,50);  ME.B = sample(1:100,50)
ME.C = sample(1:100,50);  ME.D = sample(1:100,50)  
ME1     = data.frame(ME.A, ME.B, ME.C, ME.D)
simDatA = simulateDatExpr(ME1,1000,c(0.2,0.1,0.08,0.05,0.3), signed=TRUE)
datExpr = simDatA$datExpr+5
datExpr[1:10,]  = datExpr[1:10,]+2
datExpr[41:50,] = datExpr[41:50,]-1

# Now split up the data and plot it!
subset  = c("Random Genes", "Gene.1", "Gene.234", "Gene.56", "Gene.789")
groups  = rep(c("A","A","A","B","B","B","C","C","C","C"),5)
split   = c(rep("ZZ",10), rep("YY",10), rep("XX",10), rep("WW",10), rep("VV",10))
par(mfrow = c(1,1))
results = stratifiedBarplot(datExpr, groups, split, subset)
results

# Now plot it the other way
results = stratifiedBarplot(datExpr, split, groups, subset)

}
\keyword{misc}% __ONLY ONE__ keyword per line
