\name{TrueTrait}
\alias{TrueTrait}
\title{Estimate the true trait underlying a list of surrogate markers.}
\description{
Assume an imprecisely measured trait \code{y} that is related to  the true, unobserved trait yTRUE as follows yTRUE=y+noise where noise is assumed to have mean zero and a constant variance. Assume you have 1 or more surrogate markers for yTRUE corresponding to the columns of \code{datX}. The function implements several approaches for estimating yTRUE based on the inputs \code{y} and/or \code{datX}.
}
\usage{ TrueTrait(datX, y, datXtest=NULL, 
        corFnc = "bicor", corOptions = "use = 'pairwise.complete.obs'",
        LeaveOneOut.CV=FALSE, skipMissingVariables=TRUE, 
        addLinearModel=FALSE)
}
\arguments{
  \item{datX}{ is a vector or data frame whose columns correspond to the surrogate markers (variables) for the true underlying trait. The number of rows of \code{datX} equals the number of observations, i.e. it should equal the length of \code{y} }
  
\item{y}{is a numeric vector which specifies the observed trait.
}

\item{datXtest}{can be set as a matrix or data frame of a second, independent test data set. Its columns should correspond to those of \code{datX}, i.e. the two data sets should have the same number of columns but the number or rows (test set observations) can be different.}

\item{corFnc}{Character string specifying the correlation function to be used in the calculations. 
Recomended values are the default Pearson
correlation \code{"cor"} or biweight mid-correlation \code{"bicor"}. Additional arguments to the correlation
function can be specified using \code{corOptions}.}

\item{corOptions}{Character string giving additional arguments to the function specified in \code{corFnc}. }

\item{LeaveOneOut.CV}{logical. If TRUE then leave one out cross validation estimates will be calculated for \code{y.true1} and \code{y.true2} based on \code{datX}.}

\item{skipMissingVariables}{logical. If  TRUE then variables whose values are missing for a given observation will be skipped when estimating the true trait of that particular observation. Thus, the estimate of a particular observation are determined by all the variables whose values are non-missing.}

\item{addLinearModel}{logical. If TRUE then the function also estimates the true trait based on the predictions of the linear model \code{lm(y~., data=datX)} }

}
\details{
This R function implements formulas described in Klemera and Doubal  (2006). The assumptions underlying these formulas are described in Klemera et al. But briefly,
the function provides several estimates of the true underlying trait under the following assumptions:
1) There is a true underlying trait that affects \code{y} and a list of surrogate markers corresponding to the columns of \code{datX}.
2) There is a linear relationship between the true underlying trait and \code{y} and the surrogate markers.
3)  yTRUE =y +Noise where the Noise term has a mean of zero and a fixed variance.
4) Weighted least squares estimation is used to relate the surrogate markers to the underlying trait where the weights are proportional to 1/ssq.j where ssq.j is the noise variance of the j-th marker.

Specifically,
 output \code{y.true1} corresponds to formula 31,  \code{y.true2} corresponds to formula 25, and \code{y.true3} corresponds to formula 34. 

Although the true underlying trait yTRUE is not known, one can estimate the standard deviation between the
estimate \code{y.true2} and yTRUE using formula 33. Similarly, one can estimate the SD for the estimate
\code{y.true3} using formula 42. These estimated SDs correspond to output components 2 and 3, respectively.
These SDs are valuable since they provide a sense of how accurate the measure is.  

To estimate the correlations between \code{y} and the surrogate markers, one can specify different
correlation measures. The default method is based on the Person correlation but one can also specify the
biweight midcorrelation by choosing "bicor", see help(bicor) to learn more. 

When the \code{datX} is comprised of observations measured in different strata (e.g. different batches or
independent data sets) then one can obtain stratum specific estimates by specifying the strata using the
argument \code{Strata}. In this case, the estimation focuses on one stratum at a time.  
}

\value{A list with the following components.
\item{datEstimates}{is a data frame whose columns corresponds to estimates of the true underlying trait. The number of rows equals the number of observations, i.e. the length of \code{y}.
The first column \code{y.true1} is the average value of standardized columns of \code{datX} where standardization subtracts out the intercept term and divides by the slope of the linear regression model lm(marker~y). Since this estimate ignores the fact that the surrogate markers have different correlations with \code{y}, it is typically inferior to \code{y.true2}. 
The second column \code{y.true2} equals the weighted average value of standardized columns of \code{datX}. The standardization is described in section 2.4 of Klemera et al. The weights are proportional to r^2/(1+r^2) where r denotes the correlation between the surrogate marker and \code{y}. Since this estimate does not include \code{y} as additional surrogate marker, it may be slightly inferior to \code{y.true3}. Having said this, the difference between \code{y.true2} and \code{y.true3} is often negligible. 
An additional column called \code{y.lm} is added if code{addLinearModel=TRUE}. In this case, \code{y.lm} reports the linear model predictions.
Finally, the column \code{y.true3} is very similar to \code{y.true2} but it includes \code{y} as additional surrogate marker. It is expected to be the best estimate of the underlying true trait (see Klemera et al 2006).
}

\item{datEstimatestest}{is output only if a test data set has been specified in the argument
\code{datXtest}. In this case, it contains a data frame with columns \code{ytrue1} and \code{ytrue2}. The
number of rows equals the number of test set observations, i.e the number of rows of \code{datXtest}. Since
the value of \code{y} is not known in case of a test data set, one cannot calculate \code{y.true3}. An
additional column with linear model predictions \code{y.lm} is added if code{addLinearModel=TRUE}.  }


\item{datEstimates.LeaveOneOut.CV}{is output only if the argument \code{LeaveOneOut.CV} has been set to \code{TRUE}.
In this case, it contains a data frame with leave-one-out cross validation estimates of \code{ytrue1} and \code{ytrue2}. The number of rows equals the length of \code{y}. Since the value of \code{y} is not known in case of a test data set, one cannot calculate \code{y.true3}
}



\item{SD.ytrue2}{is a scalar. This is an estimate of the standard deviation between the estimate \code{y.true2} and the true (unobserved) yTRUE. It corresponds to formula 33.}

\item{SD.ytrue3}{is a scalar. This is an estimate of the standard deviation between \code{y.true3} and the true (unobserved) yTRUE. It corresponds to formula 42.} 


\item{datVariableInfo}{is a data frame that reports information for each variable (column of \code{datX}) when it comes to the definition of \code{y.true2}. The rows correspond to the number of variables. Columns report the variable name, the center (intercept that is subtracted to scale each variable), the scale (i.e. the slope that is used in the denominator), and finally the weights used in the weighted sum of the scaled variables.}

\item{datEstimatesByStratum}{ a data frame that will only be output if \code{Strata} is different from NULL. In this case, it is has the same dimensions as \code{datEstimates} but the estimates were calculated separately for each level of \code{Strata}.}

\item{SD.ytrue2ByStratum}{ a vector of length equal to the different levels of \code{Strata}. Each component reports the estimate of \code{SD.ytrue2} for observations in the stratum specified by unique(Strata).}

\item{datVariableInfoByStratum}{ a list whose components are matrices with variable information. Each list
component reports the variable information in the stratum specified by unique(Strata). }

}

\references{
Klemera P, Doubal S (2006) A new approach to the concept and computation of biological age. Mechanisms of
Ageing and Development 127 (2006) 240-248

Choa IH, Parka KS, Limb CJ (2010)  An Empirical Comparative Study on Validation of Biological Age Estimation Algorithms with an Application of Work Ability Index. Mechanisms of Ageing and Development
Volume 131, Issue 2, February 2010, Pages 69-78 

}



\author{Steve Horvath}



\examples{
# observed trait
y=rnorm(1000,mean=50,sd=20)
# unobserved, true trait
yTRUE =y +rnorm(100,sd=10)
# now we simulate surrogate markers around the true trait
datX=simulateModule(yTRUE,nGenes=20, minCor=.4,maxCor=.9,geneMeans=rnorm(20,50,30)  )
True1=TrueTrait(datX=datX,y=y)
datTrue=True1$datEstimates
par(mfrow=c(2,2))
for (i in 1:dim(datTrue)[[2]] ){
  meanAbsDev= mean(abs(yTRUE-datTrue[,i]))
  verboseScatterplot(datTrue[,i],yTRUE,xlab=names(datTrue)[i],  
                     main=paste(i, "MeanAbsDev=", signif(meanAbsDev,3))); 
  abline(0,1)
}
#compare the estimated standard deviation of y.true2
True1[[2]]
# with the true SD
sqrt(var(yTRUE-datTrue$y.true2))
#compare the estimated standard deviation of y.true3
True1[[3]]
# with the true SD
sqrt(var(yTRUE-datTrue$y.true3))
}

\keyword{misc}


