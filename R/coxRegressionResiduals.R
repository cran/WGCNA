## The function is currently defined as

coxRegressionResiduals = function(time,event,datCovariates=NULL) 
{
if (eval(parse(text= '!require("survival")'))) 
   stop("This function requires package survival. Please install it first.");

if ( length(time) != length(event) )  { stop("Error: The length of the vector event is unequal to the length of the time vector. In R language: length(time) != length(event)")
}
if (  is.null(datCovariates) ){
coxmodel=eval(parse(text = "survival:::coxph(Surv(time, event) ~ 1 , na.action = na.exclude)"));
  }
if (  !is.null(datCovariates) ){
if ( dim(as.matrix(datCovariates))[[1]] !=length(event) ) stop("Error: the number of rows of the input matrix datCovariates is unequal to the number of observations specified in the vector event. In R language: dim(as.matrix(datCovariates))[[1]] !=length(event)")
coxmodel=eval(parse(
    text = paste("survival:::coxph(Surv(time, event) ~ . , data=datCovariates,", 
                 "na.action = na.exclude, model = TRUE)")));
} # end of if
datResiduals=data.frame(martingale=residuals(coxmodel,type="martingale"),
                        deviance=residuals(coxmodel,type="deviance"))
datResiduals
} # end of function
