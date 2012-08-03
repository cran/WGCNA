\name{metaAnalysis}
\alias{metaAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-analysis of binary and continuous variables
}
\description{
This is a meta-analysis complement to functions \code{\link{standardScreeningBinaryTrait}} and
\code{\link{standardScreeningNumericTrait}}. Given expression (or other) data from multiple independent
data sets, and the corresponding clinical traits or outcomes, the function calculates multiple screening
statistics in each data set, then calculates meta-analysis Z scores, p-values, and optionally q-values
(False Discovery Rates). Three different ways of calculating the meta-analysis Z scores are provided: the
Stouffer method, weighted Stouffer method, and using user-specified weights.
}
\usage{
metaAnalysis(multiExpr, multiTrait, 
             binary = NULL, 
             metaAnalysisWeights = NULL, 
             corFnc = cor, corOptions = list(use = "p"), 
             getQvalues = FALSE, 
             getAreaUnderROC = FALSE,
             useRankPvalue = TRUE,
             rankPvalueOptions = list(),
             setNames = NULL, 
             kruskalTest = FALSE, var.equal = FALSE, 
             metaKruskal = kruskalTest, na.action = "na.exclude")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{
   Expression data (or other data) in multi-set format (see \code{\link{checkSets}}). A vector of lists; in
each list there must be a component named \code{data} whose content
is a matrix or dataframe or array of dimension 2.
}
  \item{multiTrait}{
   Trait or ourcome data in multi-set format. Only one trait is allowed; consequesntly, the \code{data}
component of each component list can be either a vector or a data frame (matrix, array of dimension 2).
}
  \item{binary}{
   Logical: is the trait binary (\code{TRUE}) or continuous (\code{FALSE})? If not given, the decision will
be made based on the content of \code{multiTrait}.
}
  \item{metaAnalysisWeights}{
   Optional specification of set weights for meta-analysis. If given, must be a vector of non-negative
weights, one entry for each set contained in \code{multiExpr}.
}
  \item{corFnc}{
   Correlation function to be used for screening. Should be either the default \code{\link{cor}} or its
robust alternative, \code{\link{bicor}}.
}
  \item{corOptions}{
   A named list giving extra arguments to be passed to the correlation function.
}
  \item{getQvalues}{
   Logical: should q-values (FDRs) be calculated?
}
  \item{getAreaUnderROC}{
   Logical: should area under the ROC be calculated? 
   Caution, enabling the calculation will slow the function down considerably for large data sets.
}
  \item{useRankPvalue}{ Logical: should the \code{\link{rankPvalue}} function be used to obtain alternative
meta-analysis statistics?}

  \item{rankPvalueOptions}{ Additional options for function \code{\link{rankPvalue}}. These include
\code{na.last} (default \code{"keep"}), \code{ties.method} (default \code{"average"}), 
\code{calculateQvalue} (default copied from input \code{getQvalues}), 
and \code{pValueMethod} (default \code{"all"}). 
See the help file for \code{\link{rankPvalue}} for full details.}

  \item{setNames}{
   Optional specification of set names (labels). These are used to label the corresponding components of the
output. If not given, will be taken from the \code{names} attribute of \code{multiExpr}. If
\code{names(multiExpr)} is \code{NULL}, generic names of the form \code{Set_1, Set2, ...} will be used.
}
  \item{kruskalTest}{
   Logical: should the Kruskal test be performed in addition to t-test? Only applies to binary traits.
}
  \item{var.equal}{
   Logical: should the t-test assume equal variance in both groups? If \code{TRUE}, the function will warn
the user that the returned test statistics will be different from the results of the standard
\code{\link[stats]{t.test}} function.
}
  \item{metaKruskal}{
   Logical: should the meta-analysis be based on the results of Kruskal test (\code{TRUE}) or Student t-test
(\code{FALSE})?
}
  \item{na.action}{
   Specification of what should happen to missing values in \code{\link[stats]{t.test}}.
}
}
\details{
The Stouffer method of combines Z statistics by simply taking a mean of input Z statistics and multiplying
it by \code{sqrt(n)}, where \code{n} is the number of input data sets. We refer to this method as
\code{Stouffer.equalWeights}. In general, a better (i.e., more powerful) method of combining Z statistics is
to weigh them by the number of degrees of freedom (which approximately equals \code{n}). We refer to this
method as \code{weightedStouffer}. Finally, the user can also specify custom weights, for example if a data
set needs to be downweighted due to technical concerns; however, specifying own weights by hand should be
done carefully to avoid possible selection biases.
}
\value{
Data frame with the following components:
\item{ID}{ Identifier of the input genes (or other variables) }

\item{Z.equalWeights}{ Meta-analysis Z statistics obtained using Stouffer's method with equal
weights}
\item{p.equalWeights}{ p-values corresponding to \code{Z.Stouffer.equalWeights} }
\item{q.equalWeights}{ q-values corresponding to \code{p.Stouffer.equalWeights}, only present if
\code{getQvalues} is \code{TRUE}.} 

\item{Z.RootDoFWeights}{ Meta-analysis Z statistics obtained using Stouffer's method with weights given by
the square root of the number of (non-missing) samples in each data set}
\item{p.RootDoFWeights}{ p-values corresponding to \code{Z.DoFWeights} }
\item{q.RootDoFWeights}{ q-values corresponding to \code{p.DoFWeights}, only present if
\code{getQvalues} is \code{TRUE}. }

\item{Z.DoFWeights}{ Meta-analysis Z statistics obtained using Stouffer's method with weights given by
the number of (non-missing) samples in each data set}
\item{p.DoFWeights}{ p-values corresponding to \code{Z.DoFWeights} }
\item{q.DoFWeights}{ q-values corresponding to \code{p.DoFWeights}, only present if
\code{getQvalues} is \code{TRUE}. }

\item{Z.userWeights}{ Meta-analysis Z statistics
obtained using Stouffer's method with user-defined weights. Only present if input \code{metaAnalysisWeights}
are present.}
\item{p.userWeights}{ p-values corresponding to \code{Z.userWeights} }
\item{q.userWeights}{ q-values corresponding to \code{p.userWeights}, only present if
\code{getQvalues} is \code{TRUE}. }

The next set of columns is present only if input \code{useRankPvalue} is \code{TRUE} and contain the output
of the function \code{\link{rankPvalue}} with the same column weights as the above meta-analysis. Depending
on the input options \code{calculateQvalue} and \code{pValueMethod} in \code{rankPvalueOptions}, some
columns may be missing. The following columns are calculated using equal weights for each data set.

\item{pValueExtremeRank.equalWeights}{This is the minimum between pValueLowRank and
pValueHighRank, i.e. min(pValueLow, pValueHigh)}

\item{pValueLowRank.equalWeights}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the rank method.}

\item{pValueHighRank.equalWeights}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the rank method.}

\item{pValueExtremeScale.equalWeights}{This is the minimum between pValueLowScale and
pValueHighScale, i.e. min(pValueLow, pValueHigh)}

\item{pValueLowScale.equalWeights}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the Scale method.}

\item{pValueHighScale.equalWeights}{Asymptotic p-value for observing a consistently low value across
the columns of datS based on the Scale method.}

\item{qValueExtremeRank.equalWeights}{local false discovery rate (q-value) corresponding to the p-value
pValueExtremeRank}

\item{qValueLowRank.equalWeights}{local false discovery rate (q-value) corresponding to the p-value
pValueLowRank}

\item{qValueHighRank.equalWeights}{local false discovery rate (q-value) corresponding to the p-value
pValueHighRank}

\item{qValueExtremeScale.equalWeights}{local false discovery rate (q-value) corresponding to the p-value
pValueExtremeScale}

\item{qValueLowScale.equalWeights}{local false discovery rate (q-value) corresponding to the p-value
pValueLowScale}

\item{qValueHighScale.equalWeights}{local false discovery rate (q-value) corresponding to the p-value
pValueHighScale}

\item{...}{Analogous columns calculated by weighting each input set using the square root of the number of
samples, number of samples, and user weights (if given). The corresponding column names carry the suffixes
\code{RootDofWeights}, \code{DoFWeights}, \code{userWeights}.}

The following columns contain results returned by \code{\link{standardScreeningBinaryTrait}} or
\code{\link{standardScreeningNumericTrait}} (depending on whether the input trait is binary or continuous).

For binary traits, the following information is returned for each set:

  \item{corPearson.Set_1, corPearson.Set_2,...}{Pearson correlation with a binary numeric version of the
input variable. The numeric variable equals 1 for level 1 and 2 for level 2. The levels are given by
levels(factor(y)).}  

  \item{t.Student.Set_1, t.Student.Set_2, ...}{Student t-test statistic}

  \item{pvalueStudent.Set_1, pvalueStudent.Set_2, ...}{two-sided Student t-test p-value.}

  \item{qvalueStudent.Set_1, qvalueStudent.Set_2, ...}{(if input \code{qValues==TRUE})
q-value (local false discovery rate) based on the Student T-test p-value (Storey et al 2004).}

  \item{foldChange.Set_1, foldChange.Set_2, ...}{a (signed) ratio of mean values. If the mean in the first
group (corresponding to level 1) is larger than that of the second group, it equals
meanFirstGroup/meanSecondGroup. 
But if the mean of the second group is larger than that of the first group it equals
-meanSecondGroup/meanFirstGroup (notice the minus sign).}

  \item{meanFirstGroup.Set_1, meanSecondGroup.Set_2, ...}{means of columns in input \code{datExpr} across
samples in the second group.}

  \item{SE.FirstGroup.Set_1, SE.FirstGroup.Set_2, ...}{standard errors of columns in input \code{datExpr}
across samples in the 
first group.  Recall that SE(x)=sqrt(var(x)/n) where n is the number of non-missing values of x. } 

  \item{SE.SecondGroup.Set_1, SE.SecondGroup.Set_2, ...}{standard errors of columns in input \code{datExpr}
across samples in the second group.} 

  \item{areaUnderROC.Set_1, areaUnderROC.Set_2, ...}{the area under the ROC, also known as the concordance
index or C.index.  This is a measure of discriminatory power. The measure lies between 0 and 1 where 0.5
indicates no discriminatory power. 0 indicates that the "opposite" predictor has perfect discriminatory
power. To compute it we use the function \link[Hmisc]{rcorr.cens} with \code{outx=TRUE} (from Frank Harrel's
package Hmisc).}  

  \item{nPresentSamples.Set_1, nPresentSamples.Set_2, ...}{number of samples with finite measurements for
each gene.}  

If input \code{kruskalTest} is \code{TRUE}, the following columns further summarize results of 
Kruskal-Wallis test:

  \item{stat.Kruskal.Set_1, stat.Kruskal.Set_2, ...}{Kruskal-Wallis test statistic.}

  \item{stat.Kruskal.signed.Set_1, stat.Kruskal.signed.Set_2,...}{(Warning: experimental) Kruskal-Wallis
test statistic including a sign that indicates whether the average rank is higher in second group (positive)
or first group (negative). 
} 

  \item{pvaluekruskal.Set_1, pvaluekruskal.Set_2, ...}{Kruskal-Wallis test p-value.}

  \item{qkruskal.Set_1, qkruskal.Set_2, ...}{q-values corresponding to the Kruskal-Wallis test p-value (if
input \code{qValues==TRUE}).} 

  \item{Z.Set1, Z.Set2, ...}{Z statistics obtained from \code{pvalueStudent.Set1, pvalueStudent.Set2, ...}
or from \code{pvaluekruskal.Set1, pvaluekruskal.Set2, ...}, depending on input \code{metaKruskal}.}

For numeric traits, the following columns are returned:

 \item{cor.Set_1, cor.Set_2, ...}{correlations of all genes with the trait}

 \item{Z.Set1, Z.Set2, ...}{Fisher Z statistics corresponding to the correlations}

 \item{pvalueStudent.Set_1, pvalueStudent.Set_2, ...}{Student p-values of the correlations}

 \item{qvalueStudent.Set_1, qvalueStudent.Set_1, ...}{(if input \code{qValues==TRUE}) q-values of the
correlations calculated from the p-values}

 \item{AreaUnderROC.Set_1, AreaUnderROC.Set_2, ...}{area under the ROC}

 \item{nPresentSamples.Set_1, nPresentSamples.Set_2, ...}{number of samples present for the calculation of each association. }

}
\references{
For Stouffer's method, see

Stouffer, S.A., Suchman, E.A., DeVinney, L.C., Star, S.A. & Williams, R.M. Jr. 1949. The American
Soldier, Vol. 1: Adjustment during Army Life. Princeton University Press, Princeton. 

A discussion of weighted Stouffer's method can be found in 

Whitlock, M. C., Combining probability from independent tests: the weighted Z-method is superior to Fisher's
approach, Journal of Evolutionary Biology 18:5 1368 (2005)

}
\author{
Peter Langfelder
}

\seealso{
\code{\link{standardScreeningBinaryTrait}}, \code{\link{standardScreeningNumericTrait}} for screening
functions for individual data sets
}
\keyword{misc}
