\name{standardScreeningNumericTrait}
\alias{standardScreeningNumericTrait}
\title{
  Standard screening for numeric traits
}
\description{
Standard screening for numeric traits based on Pearson correlation.
}
\usage{
standardScreeningNumericTrait(datExpr, yNumeric, corFnc = cor,
                              corOptions = list(use = 'p'),
                              alternative = c("two.sided", "less", "greater"),
                              qValues = TRUE,
                              areaUnderROC = TRUE)
}
\arguments{
  \item{datExpr}{
data frame containing expression data (or more generally variables to be screened), with rows corresponding
to samples and columns to genes (variables)
}
  \item{yNumeric}{
a numeric vector giving the trait measurements for each sample
}
  \item{corFnc}{ correlation function. 
Defaults to Pearson correlation but can also be \code{\link{bicor}}. }

  \item{corOptions}{ list specifying additional arguments to be passed to the correlation function given
by \code{corFnc}. }

  \item{alternative}{alternative hypothesis for the correlation test}

  \item{qValues}{ logical: should q-values be calculated?}

  \item{areaUnderROC}{ logical: should are under the receiver-operating curve be calculated?}
}
\details{
The function calculates the correlations, associated p-values, area under the ROC, and q-values
}
\value{

Data frame with the following components:

 \item{ID }{Gene (or variable) identifiers copied from \code{colnames(datExpr)}}
 \item{cor}{correlations of all genes with the trait}
 \item{Z}{Fisher Z statistics corresponding to the correlations}
 \item{pvalueStudent }{Student p-values of the correlations}
 \item{qvalueStudent }{(if input \code{qValues==TRUE}) q-values of the correlations calculated from the p-values}
 \item{AreaUnderROC }{(if input \code{areaUnderROC==TRUE}) area under the ROC}
 \item{nPresentSamples}{number of samples present for the calculation of each association. }

}
\author{
Steve Horvath
}
\seealso{
\code{\link{standardScreeningBinaryTrait}}, \code{\link{standardScreeningCensoredTime}}
}
\keyword{misc}
