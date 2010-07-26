\name{standardScreeningNumericTrait}
\alias{standardScreeningNumericTrait}
\title{
  Standard screening for numeric traits
}
\description{
Standard screening for numeric traits based on Pearson correlation.
}
\usage{
standardScreeningNumericTrait(datExpr, yNumeric)
}
\arguments{
  \item{datExpr}{
data frame containing expression data (or more generally variables to be screened), with rows corresponding
to samples and columns to genes (variables)
}
  \item{yNumeric}{
a numeric vector giving the trait measurements for each sample
}
}
\details{
The function calculates the correlations, associated p-values, area under the ROC, and q-values
}
\value{

List with the following components:

 \item{ID }{Gene (or variable) identifiers copied from \code{colnames(datExpr)}}
 \item{corPearson }{Pearson correlations of all genes with the trait}
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
