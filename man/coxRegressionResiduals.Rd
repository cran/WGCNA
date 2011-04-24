\name{coxRegressionResiduals}
\alias{coxRegressionResiduals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Deviance- and martingale residuals from a Cox regression model
}
\description{
The function inputs a censored time variable which is specified by two input variables \code{time} and \code{event}.
It outputs i) the martingale residual and ii) deviance residual corresponding to a Cox regression model. 
By default, the Cox regression model is an intercept only Cox regression model. But optionally, the user can input covariates using the argument \code{datCovariates}.
The function makes use of the coxph function in the survival library. 
See \code{help(residuals.coxph)} to learn more.
}
\usage{
coxRegressionResiduals(time, event, datCovariates = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{ is a numeric variable that contains follow up time or time to event.
%%     ~~Describe \code{time} here~~
}
  \item{event}{ is a binary variable that takes on values 1 and 0. 1 means that the event took place (e.g. person died, or tumor recurred). 0 means censored, i.e. event has not yet been observed or loss to follow up.
}
  \item{datCovariates}{ a data frame whose columns correspond to covariates that should be used in the Cox regression model. By default, the only covariate the intercept term 1.
}
}
\details{
Residuals are often used to investigate the lack of fit of a model.  For Cox regression, there is no easy
analog to the usual "observed minus predicted" residual of linear regression. Instead, several specialized residuals have been proposed for Cox regression analysis. The function calculates residuals that are well defined for an intercept only Cox regression model: the martingale and deviance residuals (Therneau et al 1990). The martingale residual of a subject (person) specifies excess failures beyond the expected baseline hazard. 
For example, a subject who was censored at 3 years, and whose predicted cumulative hazard at 3 years was 30% has a martingale residual of 0-.30 = -0.30  
Another subject who had an event at 10 years, and whose predicted cumulative hazard at 10 years  was 60% has a martingale residual of 1-.60 = 0.40.
Since martingale residuals are not symmetrically distributed, even when the fitted model is correct, it is often advantageous to transform them into more symmetrically distributed residuals: deviance residuals. 
Thus, deviance residuals are defined as transformations of the martingale residual and the event variable. Deviance residuals are often symmetrically distributed around zero
Deviance Residuals are similar to residuals from ordinary linear regression in that they are symmetrically distributed around 0 and have standard deviation of 1.0. .  A subjects with a large deviance residual is poorly predicted by the model, i.e. is different from the baseline cumulative hazard.
A negative value indicates a longer than expected survival time.
When covariates are specified in \code{datCovariates}, then one can plot deviance (or martingale) residuals against the covariates. Unusual patterns may indicate poor fit of the Cox model.
Cryptic comments: Deviance (or martingale) residuals can sometimes be used as (uncensored) quantitative variables instead of the original time censored variable.
For example, they could be used as outcome in a regression tree or regression forest predictor.
}
\value{
It outputs a data frame with 2 columns. The first and second column correspond to martingale and deviance residuals respectively.

}
\references{
%% ~put references to the literature/web site here ~
Thereneau TM, Grambsch PM, Fleming TR (1990) Martingale-based residuals for survival models. Biometrika (1990), 77, 1, pp. 147-60

}
\author{
Steve Horvath
}
\note{
This function can be considered a wrapper of the coxph function. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
library(survival)
# simulate time and event data
time1=sample(1:100)
event1=sample(c(1,0), 100,replace=TRUE)

event1[1:5]=NA
time1[1:5]=NA
# no covariates
datResiduals= coxRegressionResiduals(time=time1,event=event1)

# now we simulate a covariate
z= rnorm(100)
cor(datResiduals,use="p")
datResiduals=coxRegressionResiduals(time=time1,event=event1,datCovariates=data.frame(z))
cor(datResiduals,use="p")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}% __ONLY ONE__ keyword per line
