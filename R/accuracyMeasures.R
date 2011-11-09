accuracyMeasures = function(tab){
if (  dim(tab)[[1]] !=2 |  dim(tab)[[2]] !=2 ) { stop("The input table is not a 2x2 table. It should be a 2x2 table, see the helpfile. ")}
if (  sum(is.na(tab) ) ) {warning("The input table is very strange. It contains missing values, i.e. some cell enties equal NA. Suggestion: check whether NA should be coded as 0.")}
is.wholenumber =function(x, tol = .Machine$double.eps^0.5) { abs(x - round(x)) < tol }
if (  sum( !is.wholenumber(tab), na.rm=T  ) >0) {warning("STRONG WARNING:
The input table contains non-integers, which does not make sense.")}

if (  sum( tab<0, na.rm=T  ) >0) {stop("The input table must be wrong since
it contains negative numbers. In R language:  sum( tab<0, na.rm=T )>0")}
num1=sum(diag(tab),na.rm=T)
denom1=sum(tab,na.rm=T)
if (denom1==0) {warning("The input table is strange. The entries sum to zero, i.e. there are zero observations. In R language: sum(tab,na.rm=T)==0.")}
TP=tab[1,1]
FP=tab[1,2]
FN=tab[2,1]
TN=tab[2,2]

error.rate= ifelse(denom1==0,NA, 1-num1/denom1)
Accuracy= ifelse(denom1==0,NA,  num1/denom1 )
Specificity= ifelse(FP + TN==0, NA,  TN / (FP + TN) )
Sensitivity= ifelse(TP + FN==0, NA,  TP / (TP + FN) )
NegativePredictiveValue= ifelse(FN + TN==0,NA,  TN / (FN + TN) )
PositivePredictiveValue=ifelse(TP + FP==0,NA,    TP / (TP + FP) )
FalsePositiveRate = 1 - Specificity 
FalseNegativeRate = 1 - Sensitivity 
Power = Sensitivity 
LikelihoodRatioPositive = ifelse(1 - Specificity==0,NA, Sensitivity / (1 - Specificity) )
LikelihoodRatioNegative = ifelse(Specificity==0, NA,  (1 - Sensitivity) / Specificity )
NaiveErrorRate = ifelse(denom1==0,NA,   min(c(tab[1,1]+ tab[2,1] , tab[1,2]+ tab[2,2] ))/denom1   )
datout=data.frame(Measure=
c("Error.Rate","Accuracy", "Specificity","Sensitivity","NegativePredictiveValue","PositivePredictiveValue","FalsePositiveRate","FalseNegativeRate","Power","LikelihoodRatioPositive","LikelihoodRatioNegative", "NaiveErrorRate"),
Value=c(error.rate,Accuracy, Specificity,Sensitivity,NegativePredictiveValue,PositivePredictiveValue,FalsePositiveRate,FalseNegativeRate,Power,LikelihoodRatioPositive,LikelihoodRatioNegative,NaiveErrorRate))
datout
  }
