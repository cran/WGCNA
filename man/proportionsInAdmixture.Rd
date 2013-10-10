\name{proportionsInAdmixture}
\alias{proportionsInAdmixture}
\title{Estimate the proportion of pure populations in an admixed population based on marker expression
values.  } 
\description{
Assume that \code{datE.Admixture} provides the expression values from a mixture of cell types (admixed
population) and you want to estimate the proportion of each pure cell type in the mixed samples (rows of
\code{datE.Admixture}).  The function allows you to do this as long as you provide a data frame
\code{MarkerMeansPure} that reports the mean expression values of markers in each of the pure cell types.  
}
\usage{
proportionsInAdmixture(
  MarkerMeansPure, 
  datE.Admixture, 
  calculateConditionNumber = FALSE, 
  coefToProportion = TRUE)
}
\arguments{
  \item{MarkerMeansPure}{ is a data frame whose first column reports the name of the marker and the
remaining columns report the mean values of the markers in each of the pure populations. The function will
estimate the proportion of pure cells which correspond to columns 2 through of
\code{dim(MarkerMeansPure)[[2]]} of \code{MarkerMeansPure}. Rows that contain missing values (NA) will be
removed.  } 
  
\item{datE.Admixture}{is a data frame of expression data, e.g. the columns of \code{datE.Admixture} could
correspond to thousands of genes.  The rows of \code{datE.Admixture} correspond to the admixed samples for
which the function estimates the proportions of pure populations.  Some of the markers specified in the
first column of \code{MarkerMeansPure} should correspond to column names of \code{datE.Admixture}. 
}
  \item{calculateConditionNumber}{logical. Default is FALSE. If set to TRUE then it uses the \code{kappa}
function to  calculates the condition number of the matrix \code{MarkerMeansPure[,-1]}.  This allows one to
determine whether the linear model for estimating the proportions is well specified. Type \code{help(kappa)}
to learn more.  \code{kappa()} computes by default (an estimate of) the 2-norm condition number of a matrix
or of the R matrix of a QR decomposition, perhaps of a linear fit. 
} 
  \item{coefToProportion}{logical. By default, it is set to TRUE. When estimating the proportions the
function fits a multivariate linear model. Ideally, the coefficients of the linear model correspond to the
proportions in the admixed samples. But sometimes the coefficients take on negative values or do not sum to
1. If \code{coefToProportion=TRUE} then negative coefficients will be set to 0 and the remaining
coefficients will be scaled so that they sum to 1. 
}
}
\details{The methods implemented in this function were motivated by
the gene expression deconvolution approach described by Abbas et al (2009), Lu  et al (2003), Wang  et al (2006). This  approach can be used to predict the proportions of (pure) cells in a complex tissue, e.g. the proportion of blood cell types in whole blood. To define the markers, you may need to have expression data from pure populations. Then you can define markers based on a significant t-test or ANOVA across the pure populations. Next use the pure population data to estimate corresponding mean expression values. Hopefully, the array platforms and normalization methods for \code{datE.MarkersAdmixtureTranspose}  and \code{MarkerMeansPure} are comparable. When dealing with Affymetrix data: we have successfully used it on untransformed MAS5 data.
For statisticians: To estimate the proportions, we use the coefficients
of a linear model. Specifically: 
\code{datCoef= t(lm(datE.MarkersAdmixtureTranspose ~MarkerMeansPure[,-1])$coefficients[-1,])}
where \code{datCoef} is a matrix whose rows correspond to the mixed samples (rows of \code{datE.Admixture})  and the columns correspond to pure populations (e.g. cell types), i.e. the columns of \code{MarkerMeansPure[,-1]}.
More details can be found in Abbas et al (2009).
}
\value{A list with the following components
\item{PredictedProportions}{data frame that contains the predicted proportions. The rows of \code{PredictedProportions} correspond to the admixed samples, i.e. the rows of \code{datE.Admixture}. The columns of \code{PredictedProportions}  correspond to the pure populations, i.e. the columns of \code{MarkerMeansPure[,-1].} }
\item{datCoef=datCoef}{data frame of numbers that is analogous to
\code{PredictedProportions}. In general, \code{datCoef} will only be different from \code{PredictedProportions} if \code{coefToProportion=TRUE}. See the description of \code{coefToProportion}
}
\item{conditionNumber}{This is the condition number resulting from the \code{kappa} function. See the description of calculateConditionNumber. }
\item{markersUsed}{vector of character strings that contains the subset of marker names (specified in the first column of \code{MarkerMeansPure}) that match column names of \code{datE.Admixture} and that contain non-missing pure mean values. }
}
\references{
Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z, Clark HF (2009) Deconvolution of Blood Microarray Data Identifies Cellular Activation Patterns in
Systemic Lupus Erythematosus. PLoS ONE 4(7): e6098. doi:10.1371/journal.pone.0006098

Lu P, Nakorchevskiy A, Marcotte EM (2003) Expression deconvolution: a
reinterpretation of DNA microarray data reveals dynamic changes in cell
populations. Proc Natl Acad Sci U S A 100: 10370-10375.

Wang M, Master SR, Chodosh LA (2006) Computational expression
deconvolution in a complex mammalian organ. BMC Bioinformatics 7: 328.
}
\author{
Steve Horvath, Chaochao Cai
}
\note{
This function can be considered a wrapper of the \code{lm} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lm}}, \code{\link{kappa}}
}
\keyword{misc}% __ONLY ONE__ keyword per line
