\name{collapseRowsUsingKME}
\alias{collapseRowsUsingKME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Selects one representative row per group based on kME }
\description{
This function selects only the most informative probe for each gene in a kME table, only keeping the probe which has the highest kME with respect to any module in the module membership matrix.  This function is a special case of the function collapseRows.
}
\usage{
collapseRowsUsingKME(MM, Gin, Pin = NULL, kMEcols = 1:dim(MM)[2])
}
\arguments{
  \item{MM}{
A module membership (kME) table with at least a subset of the columns corresponding to kME values.
}
  \item{Gin}{
Genes labels in a 1 to 1 correspondence with the rows of MM.
}
  \item{Pin}{
If NULL (default), rownames of MM are assumed to be probe IDs.  If entered, Pin must be the same length as Gin and correspond to probe IDs for MM.
}
  \item{kMEcols}{
A numeric vector showing which columns in MM correspond to kME values.  The default is all of them. 
}
}
\value{
  \item{datETcollapsed}{
A numeric matrix with the same columns as the input matrix MM, but with rows corresponding to the genes rather than the probes.
}
  \item{group2row}{
A matrix whose rows correspond to the unique gene labels and whose 2 columns report which gene label (first column called group) is represented by what probe (second column called selectedRowID)
}
  \item{selectedRow}{
A logical vector whose components are TRUE for probes selected as representatives and FALSE otherwise. It has the same length as the vector Pin.
}
}
\author{
Jeremy Miller
}

\seealso{
  \code{\link{collapseRows}}
}
\examples{
# Example: first simulate some data
set.seed(100)
ME.A = sample(1:100,50);  ME.B = sample(1:100,50)
ME.C = sample(1:100,50);  ME.D = sample(1:100,50)  
ME1     = data.frame(ME.A, ME.B, ME.C, ME.D)
simDatA = simulateDatExpr(ME1,1000,c(0.2,0.1,0.08,0.05,0.3), signed=TRUE)
simDatB = simulateDatExpr(ME1,1000,c(0.2,0.1,0.08,0.05,0.3), signed=TRUE)
Gin     = c(colnames(simDatA$datExpr),colnames(simDatB$datExpr))
Pin     = paste("Probe",1:length(Gin),sep=".")
datExpr = cbind(simDatA$datExpr, simDatB$datExpr)
MM      = corAndPvalue(datExpr,ME1)$cor

# Now run the function and see some example output
results = collapseRowsUsingKME(MM, Gin, Pin)
head(results$MMcollapsed)
head(results$group2Row)
head(results$selectedRow)

}
\keyword{misc }% __ONLY ONE__ keyword per line
