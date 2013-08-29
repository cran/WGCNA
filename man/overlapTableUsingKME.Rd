\name{overlapTableUsingKME}
\alias{overlapTableUsingKME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determines significant overlap between modules in two networks based on kME tables.  }
\description{
Takes two sets of expression data (or kME tables) as input and returns a table listing the significant overlap between each module in each data set, as well as the actual genes in common for every module pair.  Modules can be defined in several ways (generally involving kME) based on user input.
}
\usage{
overlapTableUsingKME(
   dat1, dat2, 
   colorh1, colorh2, 
   MEs1 = NULL, MEs2 = NULL, 
   name1 = "MM1", name2 = "MM2", 
   cutoffMethod = "assigned", cutoff = 0.5, 
   omitGrey = TRUE, datIsExpression = TRUE)
}
\arguments{
  \item{dat1,dat2}{
Either expression data sets (with samples as rows and genes as columns) or module membership (kME) tables (with genes as rows and modules as columns).  Function reads these inputs based on whether datIsExpression=TRUE or FALSE.  ***Be sure that these inputs include relevant row and column names, or else the function will not work properly.***
}
  \item{colorh1,colorh2}{
Color vector (module assignments) corresponding to the genes from dat1/2.  This vector must be the same length as the Gene dimension from dat1/2.
}
  \item{MEs1,MEs2}{
If entered (default=NULL), these are the module eigengenes that will be used to form the kME tables. Rows are samples and columns are module assignments.  Note that if datIsExpression=FALSE, these inputs are ignored.
}
  \item{name1,name2}{
The names of the two data sets being compared.  These names affect the output parameters.
}
  \item{cutoffMethod}{
This variable is used to determine how modules are defined in each data set.  Must be one of four options:
(1) "assigned" -> use the module assignments in colorh (default); (2) "kME" -> any gene with kME > cutoff is
in the module; (3) "numGenes" -> the top cutoff number of genes based on kME is in the module; and (4)
"pvalue" -> any gene with correlation pvalue < cutoff is in the module (this includes both positively and
negatively-correlated genes). 
}
  \item{cutoff}{
For all cutoffMethods other than "assigned", this parameter is used as the described cutoff value.
}
  \item{omitGrey}{
If TRUE the grey modules (non-module genes) for both networks are not returned.
}
  \item{datIsExpression}{
If TRUE (default), dat1/2 is assumed to be expression data.  If FALSE, dat1/2 is assumed to be a table of kME values.
}
}
\value{
  \item{PvaluesHypergeo}{
A table of p-values showing significance of module overlap based on the hypergeometric test. Note that these p-values are not corrected for multiple comparisons.
}
  \item{AllCommonGenes}{
A character vector of all genes in common between the two data sets.
}
  \item{Genes<name1/2>}{
A list of character vectors of all genes in each module in both data sets.  All genes in the MOD module in data set MM1 could be found using "<outputVariableName>$GenesMM1$MM1_MOD"
}
  \item{OverlappingGenes}{
A list of character vectors of all genes for each between-set comparison from PvaluesHypergeo.  All genes in MOD.A from MM1 that are also in MOD.B from MM2 could be found using "<outputVariableName>$OverlappingGenes$MM1_MOD.A_MM2_MOD.B"
}
}
\author{
Jeremy Miller
}
\seealso{
  \code{\link{overlapTable}}
}
\examples{
# Example: first generate simulated data.

set.seed(100)
ME.A = sample(1:100,50);  ME.B = sample(1:100,50)
ME.C = sample(1:100,50);  ME.D = sample(1:100,50) 
ME.E = sample(1:100,50);  ME.F = sample(1:100,50) 
ME.G = sample(1:100,50);  ME.H = sample(1:100,50) 
ME1     = data.frame(ME.A, ME.B, ME.C, ME.D, ME.E)
ME2     = data.frame(ME.A, ME.C, ME.D, ME.E, ME.F, ME.G, ME.H)
simDat1 = simulateDatExpr(ME1,1000,c(0.2,0.1,0.08,0.05,0.04,0.3), signed=TRUE)
simDat2 = simulateDatExpr(ME2,1000,c(0.2,0.1,0.08,0.05,0.04,0.03,0.02,0.3), 
                          signed=TRUE)

# Now run the function using assigned genes
results = overlapTableUsingKME(simDat1$datExpr, simDat2$datExpr, 
                   labels2colors(simDat1$allLabels), labels2colors(simDat2$allLabels), 
                   cutoffMethod="assigned")
results$PvaluesHypergeo

# Now run the function using a p-value cutoff, and inputting the original MEs
colnames(ME1) = standardColors(5);  colnames(ME2) = standardColors(7)
results = overlapTableUsingKME(simDat1$datExpr, simDat2$datExpr, 
                      labels2colors(simDat1$allLabels), 
                      labels2colors(simDat2$allLabels), 
                      ME1, ME2, cutoffMethod="pvalue", cutoff=0.05)
results$PvaluesHypergeo

# Check which genes are in common between the black modules from set 1 and 
# the green module from set 2
results$OverlappingGenes$MM1_green_MM2_black
}

\keyword{misc}
