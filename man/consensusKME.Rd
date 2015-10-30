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
consensusKME(
  multiExpr,
  moduleLabels, 
  multiEigengenes = NULL, 
  consensusQuantile = 0, 
  signed = TRUE,
  useModules = NULL,
  metaAnalysisWeights = NULL,
  corAndPvalueFnc = corAndPvalue, corOptions = list(), corComponent = "cor",
  getQvalues = FALSE,
  useRankPvalue = TRUE,
  rankPvalueOptions = list(calculateQvalue = getQvalues, pValueMethod = "scale"),
  setNames = NULL, 
  excludeGrey = TRUE, 
  greyLabel = if (is.numeric(moduleLabels)) 0 else "grey")
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
  \item{metaAnalysisWeights}{
    Optional specification of meta-analysis weights for each input set. If given, must be a numeric vector
of length equal the number of input data sets (i.e., \code{length(multiExpr)}). These weights will be used
in addition to constant weights and weights proportional to number of samples (observations) in each set.
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
  \item{useRankPvalue}{ Logical: should the \code{\link{rankPvalue}} function be used to obtain alternative
meta-analysis statistics?}

  \item{rankPvalueOptions}{ Additional options for function \code{\link{rankPvalue}}. These include
\code{na.last} (default \code{"keep"}), \code{ties.method} (default \code{"average"}),
\code{calculateQvalue} (default copied from input \code{getQvalues}),
and \code{pValueMethod} (default \code{"scale"}).
See the help file for \code{\link{rankPvalue}} for full details.}
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

The function \code{corAndPvalueFnc} should return a list which at the least contains (1) a matrix 
of associations of genes and eigengenes (this component should have the name given by \code{corComponent}),
and (2) a matrix of the corresponding p-values, named "p" or "p.value". Other components are optional but
for full functionality should include
(3) \code{nObs} giving the number of observations for each association (which is the number of samples less
number of missing data - this can in principle vary from association to association), and (4) \code{Z}
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

\item{weightedAverage.equalWeights.kME1, weightedAverage.equalWeights.kME2, ...}{
Average kME in each module for each gene across the
input data sets. }

\item{weightedAverage.RootDoFWeights.kME1, weightedAverage.RootDoFWeights.kME2, ...}{
Weighted average kME in each module for each gene across the
input data sets. The weight of each data set is proportional to the square root of the 
number of samples in the set. }

\item{weightedAverage.DoFWeights.kME1, weightedAverage.DoFWeights.kME2, ...}{
Weighted average kME in each module for each gene across the
input data sets. The weight of each data set is proportional to number of samples in the set. }

\item{weightedAverage.userWeights.kME1, weightedAverage.userWeights.kME2, ...}{
(Only present if input \code{metaAnalysisWeights} is non-NULL.)
Weighted average kME in each module for each gene across the
input data sets. The weight of each data set is given in \code{metaAnalysisWeights}.}

\item{meta.Z.equalWeights.kME1, meta.Z.equalWeights.kME2, ...}{Meta-analysis Z statistic for kME in each module, 
obtained by weighing the Z scores in each set equally. Only returned if the function \code{corAndPvalueFnc}
returns the Z statistics corresponding to the correlations.}

\item{meta.Z.RootDoFWeights.kME1, meta.Z.RootDoFWeights.kME2, ...}{
Meta-analysis Z statistic for kME in each module, 
obtained by weighing the Z scores in each set by the square root of the number of
samples. Only returned if the function \code{corAndPvalueFnc}
returns the Z statistics corresponding to the correlations.}

\item{meta.Z.DoFWeights.kME1, meta.Z.DoFWeights.kME2, ...}{Meta-analysis Z statistic for kME in each module, 
obtained by weighing the Z scores in each set by the number of
samples. Only returned if the function \code{corAndPvalueFnc}
returns the Z statistics corresponding to the correlations.}

\item{meta.Z.userWeights.kME1, meta.Z.userWeights.kME2, ...}{Meta-analysis Z statistic for kME in each module, 
obtained by weighing the Z scores in each set by \code{metaAnalysisWeights}. 
Only returned if \code{metaAnalysisWeights} is non-NULL and the function \code{corAndPvalueFnc}
returns the Z statistics corresponding to the correlations.}

\item{meta.p.equalWeights.kME1, meta.p.equalWeights.kME2, ...}{ 
p-values obtained from the equal-weight meta-analysis Z statistics. Only returned if the function
\code{corAndPvalueFnc} returns the Z statistics corresponding to the correlations. } 

\item{meta.p.RootDoFWeights.kME1, meta.p.RootDoFWeights.kME2, ...}{ 
p-values obtained from the meta-analysis Z statistics with weights proportional to the square root of the
number of samples. Only returned if the function
\code{corAndPvalueFnc} returns the Z statistics corresponding to the correlations. } 

\item{meta.p.DoFWeights.kME1, meta.p.DoFWeights.kME2, ...}{ 
p-values obtained from the degree-of-freedom weight meta-analysis Z statistics. Only returned if the function
\code{corAndPvalueFnc} returns the Z statistics corresponding to the correlations. } 

\item{meta.p.userWeights.kME1, meta.p.userWeights.kME2, ...}{ 
p-values obtained from the user-supplied weight meta-analysis Z statistics. Only returned if
\code{metaAnalysisWeights} is non-NULL and the function
\code{corAndPvalueFnc} returns the Z statistics corresponding to the correlations. } 

\item{meta.q.equalWeights.kME1, meta.q.equalWeights.kME2, ...}{ 
q-values obtained from the equal-weight meta-analysis p-values. Only present if
\code{getQvalues} is \code{TRUE} and the function \code{corAndPvalueFnc} 
returns the Z statistics corresponding to the kME values.}

\item{meta.q.RootDoFWeights.kME1, meta.q.RootDoFWeights.kME2, ...}{ 
q-values obtained from the meta-analysis p-values with weights proportional to the square root of the 
number of samples. Only present if
\code{getQvalues} is \code{TRUE} and the function \code{corAndPvalueFnc} 
returns the Z statistics corresponding to the kME values.}

\item{meta.q.DoFWeights.kME1, meta.q.DoFWeights.kME2, ...}{ 
q-values obtained from the degree-of-freedom weight meta-analysis p-values. Only present if
\code{getQvalues} is \code{TRUE} and the function \code{corAndPvalueFnc} 
returns the Z statistics corresponding to the kME values.}

\item{meta.q.userWeights.kME1, meta.q.userWeights.kME2, ...}{ 
q-values obtained from the user-specified weight meta-analysis p-values. Only present if
\code{metaAnalysisWeights} is non-NULL, 
\code{getQvalues} is \code{TRUE} and the function \code{corAndPvalueFnc} 
returns the Z statistics corresponding to the kME values.}

The next set of columns contain the results of function \code{\link{rankPvalue}} and are only present if
input \code{useRankPvalue} is \code{TRUE}. Some columns may be missing depending on the options specified in
\code{rankPvalueOptions}. We explicitly list columns that are based on weighing each set equally; names of
these columns carry the suffix \code{.equalWeights}

\item{pValueExtremeRank.ME1.equalWeights, pValueExtremeRank.ME2.equalWeights, ...}{
This is the minimum between pValueLowRank and
pValueHighRank, i.e. min(pValueLow, pValueHigh)}

\item{pValueLowRank.ME1.equalWeights, pValueLowRank.ME2.equalWeights, ...}{
Asymptotic p-value for observing a consistently low value
across the columns of datS based on the rank method.}

\item{pValueHighRank.ME1.equalWeights, pValueHighRank.ME2.equalWeights, ...}{
Asymptotic p-value for observing a consistently low value
across the columns of datS based on the rank method.} 

\item{pValueExtremeScale.ME1.equalWeights, pValueExtremeScale.ME2.equalWeights, ...}{
This is the minimum between pValueLowScale and
pValueHighScale, i.e. min(pValueLow, pValueHigh)}

\item{pValueLowScale.ME1.equalWeights, pValueLowScale.ME2.equalWeights, ...}{
Asymptotic p-value for observing a consistently low value
across the columns of datS based on the Scale method.} 

\item{pValueHighScale.ME1.equalWeights, pValueHighScale.ME2.equalWeights, ...}{
Asymptotic p-value for observing a consistently low
value across the columns of datS based on the Scale method.} 

\item{qValueExtremeRank.ME1.equalWeights, qValueExtremeRank.ME2.equalWeights, ...}{
local false discovery rate (q-value) corresponding
to the p-value pValueExtremeRank} 

\item{qValueLowRank.ME1.equalWeights, qValueLowRank.ME2.equalWeights, ...}{
local false discovery rate (q-value) corresponding to the
p-value pValueLowRank} 

\item{qValueHighRank.ME1.equalWeights, lueHighRank.ME2.equalWeights, ...}{
local false discovery rate (q-value) corresponding to the
p-value pValueHighRank} 

\item{qValueExtremeScale.ME1.equalWeights, qValueExtremeScale.ME2.equalWeights, ...}{
local false discovery rate (q-value)
corresponding to the p-value pValueExtremeScale}

\item{qValueLowScale.ME1.equalWeights, qValueLowScale.ME2.equalWeights, ...}{
local false discovery rate (q-value) corresponding to the
p-value pValueLowScale}

\item{qValueHighScale.ME1.equalWeights,qValueHighScale.ME2.equalWeights, ...}{
local false discovery rate (q-value) corresponding to
the p-value pValueHighScale}

\item{...}{Analogous columns corresponding to weighing individual sets by the square root of the number of
samples, by number of samples, and by user weights (if given). The corresponding column name suffixes are 
\code{.RootDoFWeights}, \code{.DoFWeights}, and \code{.userWeights}.}

The following set of columns summarize kME in individual input data sets.

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

