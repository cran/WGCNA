\name{standardScreeningBinaryTrait}
\Rdversion{1.1}
\alias{standardScreeningBinaryTrait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard screening for binatry traits
}
\description{
The function standardScreeningBinaryTrait computes widely used statistics for relating the columns of the
input data frame (argument datE) to a binary sample trait (argument y). The statistics include Student
t-test p-value and the corresponding local false discovery rate (known as q-value, Storey et al 2004),
the fold change, the area under the ROC curve (also known as C-index), mean values etc. If the input
option KruskalTest is set to TRUE, it also computes the Kruskal Wallist test p-value and corresponding
q-value. The Kruskal Wallis test is a non-parametric, rank-based group comparison test.
}
\usage{
standardScreeningBinaryTrait(datExpr, y, kruskalTest = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{
    a data frame or matrix whose columns will be related to the binary trait
}
  \item{y}{
    a binary vector whose length (number of components) equals the number of rows of datE
}
  \item{kruskalTest}{
    logical: should the Kruskal test be performed?
}
}
\value{
  A data frame whose rows correspond to the columns of datE and whose
columns report
  \item{ID}{column names of the input \code{datExpr}.}
  \item{corPearson}{pearson correlation with a binary numeric version of the input variable. The numeric
variable equals 1 for level 1 and 2 for level 2. The levels are given by levels(factor(y)).}
  \item{pvalueStudent}{two-sided Student t-test p-value.}
  \item{qvalueStudent}{q-value (local false discovery rate) based on the Student T-test p-value (Storey
et al 2004).}
  \item{foldChange}{a (signed) ratio of mean values. If the mean in the first group (corresponding to
level 1) is larger than that of the second group, it equals meanFirstGroup/meanSecondGroup.
But if the mean of the second group is larger than that of the first group it equals
-meanSecondGroup/meanFirstGroup (notice the minus sign).}
  \item{meanFirstGroup}{means of columns in input \code{datExpr} across samples in the first group.}
  \item{meanSecondGroup}{means of columns in input \code{datExpr} across samples in the second group.}
  \item{areaUnderROC}{the area under the ROC, also known as the concordance index or C.index. This is a
measure of discriminatory power. The measure lies between 0 and 1 where 0.5 indicates no discriminatory
power. 0 indicates that the "opposite" predictor has perfect discriminatory power. To compute it we use
the function rcorr.cens with outx=T (from Frank Harrel's package Hmisc).}
}
\references{
Storey JD, Taylor JE, and Siegmund D. (2004) Strong control, conservative point estimation, and
simultaneous conservative consistency of false discovery rates: A unified approach. Journal of the Royal
Statistical Society, Series B, 66: 187-205. 
}
\author{
Steve Horvath
}
\keyword{ misc }% __ONLY ONE__ keyword per line
