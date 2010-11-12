\name{standardScreeningNumericTrait}
\alias{standardScreeningNumericTrait}
\title{
  Standard screening for numeric traits
}
\description{
Standard screening for numeric traits based on Pearson correlation.
}
\usage{
standardScreeningNumericTrait(datExpr, yNumeric,
                              corFnc = "cor", corOptions = "use = 'p'")
}
\arguments{
  \item{datExpr}{
data frame containing expression data (or more generally variables to be screened), with rows corresponding
to samples and columns to genes (variables)
}
  \item{yNumeric}{
a numeric vector giving the trait measurements for each sample
}
  \item{corFnc}{ character string specifying the function to be used to calculate co-expression
similarity for correlation networks.
Defaults to Pearson correlation. Any function returning values between -1 and 1 can be used. }

  \item{corOptions}{ character string specifying additional arguments to be passed to the function given
by \code{corFnc}. Use \code{"use = 'p', method = 'Spearman'"} to obtain Spearman correlation.   }

}
\details{
The function calculates the correlations, associated p-values, area under the ROC, and q-values
}
\value{

List with the following components:

 \item{ID }{Gene (or variable) identifiers copied from \code{colnames(datExpr)}}
 \item{cor}{correlations of all genes with the trait}
 \item{pvalueStudent }{Student p-values of the correlations}
 \item{qvalueStudent }{q-values of the correlations calculated from the p-values}
 \item{AreaUnderROC }{area under the ROC}

}
\author{
Steve Horvath
}
\seealso{
\code{\link{standardScreeningBinaryTrait}}, \code{\link{standardScreeningCensoredTime}}
}
\keyword{misc}
