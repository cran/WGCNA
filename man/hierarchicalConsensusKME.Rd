\name{hierarchicalConsensusKME}
\alias{hierarchicalConsensusKME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of measures of fuzzy module membership (KME) in hierarchical consensus modules
}
\description{
This function calculates several measures of fuzzy module membership in hiearchical consensus modules.
}
\usage{
hierarchicalConsensusKME(
   multiExpr,
   moduleLabels,
   multiEigengenes = NULL,
   consensusTree,
   signed = TRUE,
   useModules = NULL,
   metaAnalysisWeights = NULL,
   corAndPvalueFnc = corAndPvalue, corOptions = list(),
   corComponent = "cor", getFDR = FALSE,
   useRankPvalue = TRUE,
   rankPvalueOptions = list(calculateQvalue = getFDR, pValueMethod = "scale"),
   setNames = NULL, excludeGrey = TRUE,
   greyLabel = if (is.numeric(moduleLabels)) 0 else "grey",
   reportWeightType = NULL,
   getOwnModuleZ = TRUE,
   getBestModuleZ = TRUE,
   getOwnConsensusKME = TRUE,
   getBestConsensusKME = TRUE,
   getAverageKME = FALSE,
   getConsensusKME = TRUE,

   getMetaP = FALSE,
   getMetaFDR = getMetaP && getFDR,

   getSetKME = TRUE,
   getSetZ = FALSE,
   getSetP = FALSE,
   getSetFDR = getSetP && getFDR,

   includeID = TRUE,
   additionalGeneInfo = NULL,
   includeWeightTypeInColnames = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{multiExpr}{ Expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }
    \item{moduleLabels}{
A vector with one entry per column (gene or probe) in \code{multiExpr}, giving the module labels. 
}
  \item{multiEigengenes}{
Optional specification of module eigengenes of the modules (\code{moduleLabels}) in data sets within
\code{multiExpr}. If not given, will be calculated.
}
  \item{consensusTree}{
A list specifying the consensus calculation. See details.
}
  \item{signed}{
Logical: should module membership be considered singed? Signed membership should be used for signed (including
signed hybrid) networks and means that negative module membership means the
gene is not a member of the module.
    In other words, in signed networks 
    negative kME values are not considered significant and the corresponding p-values will be one-sided. In
    unsigned networks, negative kME values are considered significant and the corresponding
    p-values will be two-sided.
}
  \item{useModules}{
Optional vector specifying which modules should be used. Defaults to all modules except the unassigned module.
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
  \item{getFDR}{
Logical: should FDR be calculated?
}
  \item{useRankPvalue}{ Logical: should the \code{\link{rankPvalue}} function be used to obtain alternative
meta-analysis statistics?}
  \item{rankPvalueOptions}{ Additional options for function \code{\link{rankPvalue}}. These include
\code{na.last} (default \code{"keep"}), \code{ties.method} (default \code{"average"}),
\code{calculateQvalue} (default copied from input \code{getQvalues}),
and \code{pValueMethod} (default \code{"scale"}).
See the help file for \code{\link{rankPvalue}} for full details.}
  \item{setNames}{
   Names for the input sets. If not given, will be taken from \code{names(multiExpr)}. If those are
\code{NULL} as well, the names will be \code{"Set_1", "Set_2", ...}.
}
   \item{excludeGrey}{
   logical: should the grey module be excluded from the kME tables? Since the grey module is typically not a
real module, it makes little sense to report kME values for it.
}
   \item{greyLabel}{
   label that labels the grey module.
}
  \item{reportWeightType}{
One of \code{"equal", "rootDoF", "DoF", "user"}. Indicates which of the weights should be reported in the
output. If not given, all available weight types will be reported; this always includes \code{"equal",
"rootDoF", "DoF"}, while \code{"user"} weights are reported if \code{metaAnalysisWeights} above is given.
}
  \item{getOwnModuleZ}{
Logical: should meta-analysis Z statistic in own module be returned as a column of the output?
}
  \item{getBestModuleZ}{
Logical: should highest meta-analysis Z statistic across all modules and the corresponding module be returned
as columns of the output?
}
  \item{getOwnConsensusKME}{
Logical: should consensus KME (eigengene-based connectivity) statistic in own module be returned as a 
column of the output?
}
  \item{getBestConsensusKME}{
Logical: should highest consensus KME across all modules and the corresponding module be returned 
as columns of the output?
}
  \item{getAverageKME}{
Logical: Should average KME be calculated? 
}
  \item{getConsensusKME}{
Logical: should consensus KME be calculated?
}
  \item{getMetaP}{
Logical: should meta-analysis p-values corresponding to the KME meta-analysis Z statistics be calculated?
}
  \item{getMetaFDR}{
Logical: should FDR estimates for the meta-analysis p-values corresponding to the KME meta-analysis Z
statistics be calculated? 
}
  \item{getSetKME}{
Logical: should KME values for individual sets be returned?
}
  \item{getSetZ}{
Logical: should Z statistics corresponding to KME for individual sets be returned?
}
  \item{getSetP}{
Logical: should p values corresponding to KME for individual sets be returned?
}
  \item{getSetFDR}{
Logical: should FDR estimates corresponding to KME for individual sets be returned?
}
  \item{includeID}{
Logical: should gene ID (taken from column names of \code{multiExpr}) be included as the first column in
the output?
}
  \item{additionalGeneInfo}{
Optional data frame with rows corresponding to genes in \code{multiExpr} that should be included as part of
the output.
}
  \item{includeWeightTypeInColnames}{
Logical: should weight type (\code{"equal", "rootDoF", "DoF", "user"}) be included in appropriate
meta-analysis column names?
}
}
\details{
This function calculates several measures of (hierarchical) consensus KME (eigengene-based intramodular
connectivity or fuzzy module membership) for all genes in all modules.

First, it calculates the meta-analysis Z statistics for correlations
between genes and module eigengenes; this is known as the consensus module membership Z statistic. The
meta-analysis weights can be specified by the user either explicitly or implicitly ("equal", "RootDoF" or
"DoF"). 

Second, it can calculate the consensus KME, i.e., the hierarchical consensus of the KMEs (correlations with
eigengenes) across the individual sets. The consensus calculation is specified in the argument
\code{consensusTree}; 
typically, the \code{consensusTree} used here will be the same as the one used for the actual consensus
network construction and module identification.  
See \code{\link{newConsensusTree}} for details on how to specify consensus trees.

Third, the function can also calculate the (weighted) average KME using the meta-analysis weights; the average
KME can be interpreted as the meta-analysis of the KMEs in the individual sets. This is related to but
somewhat distinct from the meta-analysis Z statistics.

In addition to these, optional output also includes, for each gene, KME values in the module to which the gene
is assigned as well as the maximum KME values and modules for
which the maxima are attained. For most genes, the assigned module will be the one with highest KME values,
but for some genes the assigned module and module of maximum KME may be different.

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
Data frame with the following components, some of which may be missing depending on input options (for easier
readability the order here is not the same as in the actual output): 

\item{ID}{Gene ID, taken from the column names of the first input data set}

If given, a copy of \code{additionalGeneInfo}.

\item{Z.kME.inOwnModule}{Meta-analysis Z statistic for membership in assigned module.}

\item{maxZ.kME}{Maximum meta-analysis Z statistic for membership across all modules.}
\item{moduleOfMaxZ.kME}{Module in which the maximum meta-analysis Z statistic is attained. }

\item{consKME.inOwnModule}{Consensus KME in assigned module.}

\item{maxConsKME}{Maximum consensus KME across all modules.}
\item{moduleOfMaxConsKME}{Module in which the maximum consensus KME is attained.}


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
based on the rank method.}

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
\author{
Peter Langfelder
}

\seealso{
\code{\link{signedKME}} for eigengene based connectivity in a single data set.
\code{\link{corAndPvalue}}, \code{\link{bicorAndPvalue}} for two alternatives for calculating correlations and the 
corresponding p-values and Z scores. Both can be used with this function.
\code{\link{newConsensusTree}} for more details on hierarchical consensus trees and calculations.
}
\keyword{misc}
