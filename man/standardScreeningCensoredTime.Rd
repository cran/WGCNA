\name{standardScreeningCensoredTime}
\Rdversion{1.1}
\alias{standardScreeningCensoredTime}
\title{
Standard Screening with regard to a Censored Time Variable 
}
\description{
  The function standardScreeningCensoredTime computes association measures between the columns of the input
data datE and a censored time variable (e.g. survival time). The censored time is specified using two input
variables "time" and "event". The event variable is binary where 1 indicates that the event took place (e.g.
the person died) and 0 indicates censored (i.e. lost to follow up).  The function fits univariate Cox
regression models (one for each column of datE) and outputs a Wald test p-value, a logrank p-value,
corresponding local false discovery rates (known as q-values, Storey et al 2004), hazard ratios. Further it
reports the concordance index (also know as area under the ROC curve) and optionally results from
dichotomizing the columns of datE.

}
\usage{
standardScreeningCensoredTime(
   time, 
   event, 
   datExpr, 
   percentiles = seq(from = 0.1, to = 0.9, by = 0.2), 
   dichotomizationResults = FALSE, 
   qValues = TRUE,
   fastCalculation = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
numeric variable showing time to event or time to last follow up.
}
  \item{event}{
Input variable \code{time} specifies the time to event or time to last follow up. Input variable
\code{event} indicates whether the event happend (=1) or whether there was censoring (=0). 
}
  \item{datExpr}{
a data frame or matrix whose columns will be related to the censored time.
}
  \item{percentiles}{
numeric vector which is only used when dichotomizationResults=T. Each value should lie between 0 and 1. For
each value specified in the vector percentiles, a binary vector will be defined by dichotomizing the column
value according to the corresponding quantile. Next a corresponding p-value will be calculated.
}
  \item{dichotomizationResults}{
logical. If this option is set to TRUE then the values of the columns of datE will be dichotomized and
corresponding Cox regression p-values will be calculated. 
}
  \item{qValues}{
 logical. If this option is set to TRUE (default) then q-values will be calculated for the Cox regression
p-values.
}
  \item{fastCalculation}{
logical. If set to TRUE, the function outputs correlation test p-values (and q-values) for correlating the
columns of datE with the expected hazard (if no covariate is fit). Specifically, the expected hazard is
defined as the deviance residual of an intercept only Cox regression model. The results are very similar to
those resulting from a univariate Cox model where the censored time is regressed on the columns of dat.
Specifically, this computational speed up is facilitated by the insight that the p-values resulting from a
univariate Cox regression coxph(Surv(time,event)~datE[,i]) are very similar to those from
corPvalueFisher(cor(devianceResidual,datE[,i]), nSamples).
}
}
\details{

If input option fastCalculation=TRUE, then the function outputs correlation test p-values (and q-values) for
correlating the columns of datE with the expected hazard (if no covariate is fit). Specifically, the
expected hazard is defined as the deviance residual of an intercept only Cox regression model. The results
are very similar to those resulting from a univariate Cox model where the censored time is regressed on the
columns of dat. Specifically, this computational speed up is facilitated by the insight that the p-values
resulting from a univariate Cox regression coxph(Surv(time,event)~datE[,i]) are very similar to those from
corPvalueFisher(cor(devianceResidual,datE[,i]), nSamples)


}
\value{
  If \code{fastCalculation} is \code{FALSE}, 
the function outputs a data frame whose rows correspond to the columns of datE and whose columns report 
  \item{ID}{column names of the input data datExpr.}
  \item{pvalueWald}{Wald test p-value from fitting a univariate Cox regression model where the censored time
is regressed on each column of datExpr.}
  \item{qValueWald}{local false discovery rate (q-value) corresponding to the Wald test p-value. }
  \item{pvalueLogrank}{Logrank p-value resulting from the Cox regression model. Also known as score test
p-value. For large sample sizes this sould be similar to the Wald test p-value. }
  \item{qValueLogrank}{local false discovery rate (q-value) corresponding to the Logrank test p-value. }
  \item{HazardRatio}{hazard ratio resulting from the Cox model. If the value is larger than 1, then high
values of the column are associated with shorter time, e.g. increased hazard of death. A hazard ratio equal
to 1 means no relationship between the column and time. HR<1 means that high values are associated with
longer time, i.e. lower hazard.}
  \item{CI.LowerLimitHR}{Lower bound of the 95 percent confidence interval of the hazard ratio. }
  \item{CI.UpperLimitHR}{Upper bound of the 95 percent confidence interval of the hazard ratio. }
  \item{C.index}{concordance index, also known as C-index or area under the ROC curve. Calculated with the
rcorr.cens option outx=TRUE (ties are ignored).}
  \item{MinimumDichotPvalue}{This is the smallest p-value from the dichotomization results. To see which
dichotomized variable (and percentile) corresponds to the minimum, study the following columns. }
  \item{pValueDichot0.1}{This columns report the p-value when the column is dichotomized according to the
specified percentile (here 0.1). The percentiles are specified in the input option percentiles. }
  \item{pvalueDeviance}{The p-value resulting from using a correlation test to relate the expected hazard
(deviance residual) with each (undichotomized) column of datE. Specifically, the Fisher transformation is
used to calculate the p-value for the Pearson correlation. The resulting p-value should be very similar to
that of a univariate Cox regression model.}
  \item{qvalueDeviance}{Local false discovery rate (q-value) corresponding to pvalueDeviance.}
  \item{corDeviance}{Pearson correlation between the expected hazard (deviance residual) with each
(undichotomized) column of datExpr.}
}
\author{
Steve Horvath
}
\keyword{ misc }
