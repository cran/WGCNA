\name{kMEcomparisonScatterplot}
\alias{kMEcomparisonScatterplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to plot kME values between two comparable data sets.  }
\description{
Plots the kME values of genes in two groups of expression data for each module in an inputted color vector.
}
\usage{
kMEcomparisonScatterplot(
   datExpr1, datExpr2, colorh, 
   inA = NULL, inB = NULL, MEsA = NULL, MEsB = NULL, 
   nameA = "A", nameB = "B", 
   plotAll = FALSE, noGrey = TRUE, maxPlot = 1000, pch = 19, 
   fileName = if (plotAll) paste("kME_correlations_between_",nameA,"_and_",
                                 nameB,"_all.pdf",sep="") else
                           paste("kME_correlations_between_",nameA,"_and_",
                                 nameB,"_inMod.pdf",sep=""), ...)
}
\arguments{
  \item{datExpr1}{
The first expression matrix (samples=rows, genes=columns).  This can either include only the data for group A (in which case dataExpr2 must be entered), or can contain all of the data for groups A and B (in which case inA and inB must be entered).
}
  \item{datExpr2}{
The second expression matrix, or set to NULL if all data is from same expression matrix.  If entered, datExpr2 must contain the same genes as datExpr1 in the same order.
}
  \item{colorh}{
The common color vector (module labels) corresponding to both sets of expression data.
}
  \item{inA, inB}{
Vectors of TRUE/FALSE indicating whether a sample is in group A/B, or a vector of numeric indices indicating which samples are in group A/B. If datExpr2 is entered, these inputs are ignored (thus default = NULL).  For these and all other A/B inputs, "A" corresponds to datExpr1 and "B" corresponds to datExpr2 if datExpr2 is entered; otherwise "A" corresponds to datExpr1[inA,] while "B" corresponds to datExpr1[inB,].
}
  \item{MEsA, MEsB}{
Either the module eigengenes or NULL (default) in which case the module eigengenes will be calculated.  In inputted, MEs MUST be calculated using "moduleEigengenes(<parameters>)$eigengenes" for function to work properly.
}
  \item{nameA, nameB}{
The names of these groups (defaults = "A" and "B").  The resulting file name (see below) and x and y axis labels for each scatter plot depend on these names.
}
  \item{plotAll}{
If TRUE, plot gene-ME correlations for all genes.  If FALSE, plot correlations for only genes in the plotted module (default).  Note that the output file name will be different depending on this parameter, so both can be run without overwriting results.
}
  \item{noGrey}{
If TRUE (default), the grey module genes are ignored.  This parameter is only used if MEsA and MEsB are calculated.
}
  \item{maxPlot}{
The maximum number of random genes to include (default=1000).  Smaller values lead to smaller and less cluttered plots, usually without significantly affecting the resulting correlations. This parameter is only used if plotAll=TRUE.
}
  \item{pch}{
See help file for "points". Setting pch=19 (default) produces solid circles.
}
  \item{fileName}{ Name of the file to hold the plots. Since the output format is pdf, the extension should
be .pdf .}

  \item{...}{
Other plotting parameters that are allowable inputs to verboseScatterplot.
}
}
\value{
The default output is a file called "kME_correlations_between_[nameA]_and_[nameB]_[all/inMod].pdf", where
[nameA] and [nameB] correspond to the nameA and nameB input parameters, and [all/inMod] depends on whether
plotAll=TRUE or FALSE. This output file contains all of the plots as separate pdf images, and will be
located in the current working directory.  
}
\author{
Jeremy Miller
}
\note{
The function "pdf", which can be found in the grDevices library, is required to run this function.
}
\examples{
# Example output file ("kME_correlations_between_A_and_B_inMod.pdf") using simulated data.
\dontrun{
set.seed = 100
ME=matrix(0,50,5)
for (i in 1:5) ME[,i]=sample(1:100,50)
simData1 = simulateDatExpr5Modules(MEturquoise=ME[,1],MEblue=ME[,2],
                          MEbrown=ME[,3],MEyellow=ME[,4], MEgreen=ME[,5])
simData2 = simulateDatExpr5Modules(MEturquoise=ME[,1],MEblue=ME[,2],
                          MEbrown=ME[,3],MEyellow=ME[,4], MEgreen=ME[,5])
kMEcomparisonScatterplot(simData1$datExpr,simData2$datExpr,simData1$truemodule)
}

}
\keyword{misc}
