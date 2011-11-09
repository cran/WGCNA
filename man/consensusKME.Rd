\name{consensusKME}
\alias{consensusKME}
\title{
Calculate consensus kME (eigengene-based connectivities) across multiple data sets.
}
\description{
Calculate consensus kME (eigengene-based connectivities) across multiple data sets, typically following a
consensus module analysis. 
}
\usage{
consensusKME(multiExpr, moduleLabels, 
             multiEigengenes = NULL, 
             consensusQuantile = 0, 
             signed = TRUE,
             useModules = NULL,
             countWeightPower = 1,
             corAndPvalueFnc = corAndPvalue, corOptions = list(), corComponent = "cor",
             getQvalues = FALSE,
             setNames = NULL, 
             excludeGrey = TRUE, greyLabel = ifelse(is.numeric(moduleLabels), 0, "grey"))
}
\arguments{
  \item{multiExpr}{
  Expression (or other numeric) data in a multi-set format. A vector of lists; in each list there must be a
component named `data' whose content is a matrix or dataframe or array of dimension 2.
}
  \item{moduleLabels}{
  Module labels: one label for each gene in \code{multiExpr}.
}
  \item{multiEigengenes}{
  Optional eigengenes of modules specified in \code{moduleLabels}. If not given, will be calculated from
\code{multiExpr}. 
}
  \item{signed}{
    logical: should the network be considered signed? In signed networks (\code{TRUE}), 
    negative kME values are not considered significant and the corresponding p-values will be one-sided. In
    unsigned networks (\code{FALSE}), negative kME values are considered significant and the corresponding
    p-values will be two-sided.
}
  \item{useModules}{
  Optional specification of module labels to which the analysis should be restricted. This could be useful
if there are many modules, most of which are not interesting. Note that the "grey" module cannot be used
with \code{useModules}.}

  \item{consensusQuantile}{
  Quantile for the consensus calculation. Should be a number between 0 (minimum) and 1.
}
  \item{countWeightPower}{
  This argument controls how the individual kME values are weighted in the weighted average output 
  (see Value below). The weight of each input set is proportional to number of samples in the set, 
  raised to the power \code{countWeightPower}. 
}
  \item{corAndPvalueFnc}{
  Function that calculates associations between expression profiles and eigengenes. See details.
}
  \item{corOptions}{
  List giving additional arguments to function \code{corAndPvalueFnc}. See details.
}
  \item{corComponent}{
  Name of the component of output of \code{corAndPvalueFnc} that contains the actual correlation.
}
  \item{getQvalues}{
   logical: should q-values (estimates of FDR) be calculated?
}
  \item{setNames}{
   names for the input sets. If not given, will be taken from \code{names(multiExpr)}. If those are
\code{NULL} as well, the names will be \code{"Set_1", "Set_2", ...}.
}
   \item{excludeGrey}{
   logical: should the grey module be excluded from the kME tables? Since the grey module is typically not a
real module, it makes little sense to report kME values for it.
}
   \item{greyLabel}{
   label that labels the grey module.
}
}
\details{
The function \code{corAndPvalueFnc} is currently
is expected to accept arguments \code{x} (gene expression profiles), \code{y} (eigengene expression
profiles), and \code{alternative} with possibilities at least \code{"greater", "two.sided"}. 
Any additional arguments can be passed via \code{corOptions}. 

The function \code{corAndPvalueFnc} should return a list which at the least contains a matrix 
of associations of genes and eigengenes.
This component should have the name given by \code{corComponent} below. Other components are optional but
for full functionality should include
(1) \code{nObs} giving the number of observations for each association (which is the number of samples less
number of missing data - this can in principle vary from association to association), and (2) \code{Z}
giving a Z static for each observation. If these are missing, \code{nObs} is calculated in the main
function, and calculations using the Z statistic are skipped.
}
\value{
Data frame with the following components (for easier readability the order here is not the same as in the
actual output):
\item{ID}{Gene ID, taken from the column names of the first input data set}

\item{consensus.kME.1, consensus.kME.2, ...}{Consensus kME (that is, the requested quantile of the kMEs in the
individual data sets)in each module for each gene across the input data
sets. The module labels (here 1, 2, etc.) correspond to those in \code{moduleLabels}.}

\item{average.kME.1, average.kME.2, ...}{Average kME in each module for each gene across the input data 
sets. }

\item{weightedAverage.kME.1, weightedAverage.kME.2, ...}{Average kME in each module for each gene across the
input data sets. The weight of each data set is proportional to number of samples in the set, raised to the
power \code{countWeightPower}. } 

\item{meta.Z.kME1, meta.Z.kME2, ...}{Meta-analysis Z statistic for kME in each module, 
obtained by weighing the Z scores in each set by number of
samples raised to the power \code{countWeightPower}. Only returned if the function \code{corAndPvalueFnc}
returns the Z statistics corresponding to the correlations.}

\item{meta.p.kME1, meta.p.kME2, ...}{ p-values obtained from the meta-analysis Z statistics.Only returned if
the function \code{corAndPvalueFnc} 
returns the Z statistics corresponding to the correlations. }

\item{meta.q.kME1, meta.q.kME2, ...}{ q-values obtained from the meta-analysis p-values. Only present if
\code{getQvalues} is \code{TRUE} and the function \code{corAndPvalueFnc} 
returns the Z statistics corresponding to the kME values.}

\item{kME1.Set_1, kME1.Set_2, ..., kME2.Set_1, kME2.Set_2, ...}{ kME values for each gene in each module in
each given data set. }

\item{p.kME1.Set_1, p.kME1.Set_2, ..., p.kME2.Set_1, p.kME2.Set_2, ...}{ p-values corresponding to 
kME values for each gene in each module in each given data set. }

\item{q.kME1.Set_1, q.kME1.Set_2, ..., q.kME2.Set_1, q.kME2.Set_2, ...}{ q-values corresponding to 
kME values for each gene in each module in each given data set. Only returned if \code{getQvalues} is
\code{TRUE}. }

\item{Z.kME1.Set_1, Z.kME1.Set_2, ..., Z.kME2.Set_1, Z.kME2.Set_2, ...}{ Z statistics corresponding to
kME values for each gene in each module in each given data set. Only present if the function
\code{corAndPvalueFnc}                 
returns the Z statistics corresponding to the kME values. }

}
\references{
Langfelder P, Horvath S., WGCNA: an R package for weighted correlation network analysis.
BMC Bioinformatics. 2008 Dec 29; 9:559.

}
\author{
Peter Langfelder
}

\seealso{
\link{signedKME} for eigengene based connectivity in a single data set.
\link{corAndPvalue}, \link{bicorAndPvalue} for two alternatives for calculating correlations and the
corresponding p-values and Z scores. Both can be used with this function.
}
\keyword{misc}

