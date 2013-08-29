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
standardScreeningBinaryTrait(
     datExpr, y, 
     corFnc = cor, corOptions = list(use = 'p'),
     kruskalTest = FALSE, qValues = FALSE,
     var.equal=FALSE, na.action="na.exclude",
     getAreaUnderROC = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{
    a data frame or matrix whose columns will be related to the binary trait
}
  \item{y}{
    a binary vector whose length (number of components) equals the number of rows of datE
}
  \item{corFnc}{ correlation function. Defaults to Pearson correlation. 
}
  \item{corOptions}{ a list specifying options to corFnc. An empty list must be specified as \code{list()}
(supplying \code{NULL} instead will trigger an error). 
}

  \item{kruskalTest}{
    logical: should the Kruskal test be performed?
}
  \item{qValues}{
    logical: should the q-values be calculated?
}
  \item{var.equal}{
    logical input parameter for the Student t-test. It indicates whether to treat the two variances
(corresponding to the binary grouping) are being equal. If TRUE then the pooled variance is used to estimate
the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
Warning: here the default value is TRUE which is different from the default value of t.test. Type help(t.test)
for more details. }
\item{na.action}{
    character string for the Student t-test: indicates what should happen when the data contain missing
values NAs. }
\item{getAreaUnderROC}{logical: should area under the ROC curve be calculated? The calculation slows the
function down somewhat. }

}

\value{
  A data frame whose rows correspond to the columns of datE and whose
columns report
  \item{ID}{column names of the input \code{datExpr}.}
  \item{corPearson}{pearson correlation with a binary numeric version of the input variable. The numeric
variable equals 1 for level 1 and 2 for level 2. The levels are given by levels(factor(y)).}
  \item{t.Student}{Student's t-test statistic}
  \item{pvalueStudent}{two-sided Student t-test p-value.}
  \item{qvalueStudent}{(if input \code{qValues==TRUE}) 
q-value (local false discovery rate) based on the Student T-test p-value (Storey et al 2004).}
  \item{foldChange}{a (signed) ratio of mean values. If the mean in the first group (corresponding to
level 1) is larger than that of the second group, it equals meanFirstGroup/meanSecondGroup.
But if the mean of the second group is larger than that of the first group it equals
-meanSecondGroup/meanFirstGroup (notice the minus sign).}
  \item{meanFirstGroup}{means of columns in input \code{datExpr} across samples in the first group.}
  \item{meanSecondGroup}{means of columns in input \code{datExpr} across samples in the second group.}

  \item{SE.FirstGroup}{standard errors of columns in input \code{datExpr} across samples in the first group. Recall that SE(x)=sqrt(var(x)/n) where n is the number of non-missing values of x. }
  \item{SE.SecondGroup}{standard errors of columns in input \code{datExpr} across samples in the second group.}
  \item{areaUnderROC}{the area under the ROC, also known as the concordance index or C.index. This is a
measure of discriminatory power. The measure lies between 0 and 1 where 0.5 indicates no discriminatory
power. 0 indicates that the "opposite" predictor has perfect discriminatory power. To compute it we use
the function \link[Hmisc]{rcorr.cens} with \code{outx=TRUE} (from Frank Harrel's package Hmisc).
Only present if input \code{getAreUnderROC} is \code{TRUE}.}
  \item{nPresentSamples}{number of samples with finite measurements for each gene.}

If input \code{kruskalTest} is \code{TRUE}, the following columns further summarize results of
Kruskal-Wallis test:
  \item{stat.Kruskal}{Kruskal-Wallis test statistic.}
  \item{stat.Kruskal.signed}{(Warning: experimental) Kruskal-Wallis test statistic including a sign that
indicates whether the average rank is higher in second group (positive) or first group (negative). }
  \item{pvaluekruskal}{Kruskal-Wallis test p-values.}
  \item{qkruskal}{q-values corresponding to the Kruskal-Wallis test p-value (if input \code{qValues==TRUE}).}
}

\references{
Storey JD, Taylor JE, and Siegmund D. (2004) Strong control, conservative point estimation, and
simultaneous conservative consistency of false discovery rates: A unified approach. Journal of the Royal
Statistical Society, Series B, 66: 187-205. 
}
\author{
Steve Horvath
}

\examples{

require(survival) # For is.Surv in rcorr.cens
m=50
y=sample(c(1,2),m,replace=TRUE)
datExprSignal=simulateModule(scale(y),30)
datExprNoise=simulateModule(rnorm(m),150)
datExpr=data.frame(datExprSignal,datExprNoise)

Result1=standardScreeningBinaryTrait(datExpr,y)
Result1[1:5,]



# use unequal variances and calculate q-values
Result2=standardScreeningBinaryTrait(datExpr,y, var.equal=FALSE,qValue=TRUE)
Result2[1:5,]

# calculate Kruskal Wallis test and q-values
Result3=standardScreeningBinaryTrait(datExpr,y,kruskalTest=TRUE,qValue=TRUE)
Result3[1:5,]

}

\keyword{ misc }% __ONLY ONE__ keyword per line

