
TrueTrait=function(datX, y,datXtest=NULL, corFnc = "bicor", corOptions = "use = 'pairwise.complete.obs'", LeaveOneOut.CV=FALSE,skipMissingVariables=TRUE,addLinearModel=FALSE){
datX=as.matrix(datX)
no.variables=dim(as.matrix(datX))[[2]]
datVariableInfo=data.frame(matrix(NA, nrow=no.variables, ncol=4))
names(datVariableInfo)=c("Variable","center", "scale", "weights.y.true2")
if ( is.null(colnames(datX))){datVariableInfo$Variable=1:no.variables} else { datVariableInfo$Variable=colnames(datX)}
no.observations=dim(as.matrix(datX))[[1]]
if (no.observations !=  length(y) ) {stop("The number of rows of datX does not correspond to the length of y. Consider transposing your input matrix or use a different input.") }
if (no.observations==1 ) {
warning("Only 1 observations, i.e. the length of y is 1. The function cannot be used. For your convenience, the estimated true values will be set to the input value of y.")
y.true1=y; y.true2=y;  y.true3=y}
if (no.observations>1 ) {
y.true1=rep(NA, length(y) )
y.true2=rep(NA, length(y) )
y.true3=rep(NA, length(y) )
y.lm=rep(NA, length(y) )
restNonMissingY= !is.na(y)
r.characteristic=NA
SD.ytrue2=NA
SD.ytrue3=NA
SsquaredBE=NA

if (sum(restNonMissingY,na.rm=TRUE) >3 ){
corX = parse(text = paste(corFnc, "(datX,y ",prepComma(corOptions), ")"))
rVector= as.numeric(eval(corX))
datCoef=t(coef(lm(datX~ y,na.action="na.exclude")))
datVariableInfo$center=datCoef[,1]	#intercept
datVariableInfo$scale=datCoef[,2]	#slope
datXscaled=scale(datX,center=datCoef[,1],scale=datCoef[,2] )
weights0=rVector^2/((1-rVector^2)*var(y,na.rm=TRUE)) # Steve, this is where I made the one change
weights=weights0/sum(weights0)
datVariableInfo$weights.y.true2=weights
y.true1=as.numeric(apply(as.matrix(datXscaled),1,mean))
y.true2=as.numeric(as.matrix(datXscaled)%*%weights)

if (skipMissingVariables ) {
y.true1= as.numeric(apply(as.matrix(datXscaled),1,mean,na.rm=TRUE))
weightsMatrix=matrix(weights,byrow=TRUE,nrow=dim(as.matrix(datXscaled))[[1]],ncol=length(weights) )
weightsMatrix[is.na(datXscaled)]=0
rowsum.weightsMatrix=apply(as.matrix(weightsMatrix),1,sum)
weightsMatrix=t(scale(t(as.matrix(weightsMatrix)),center=F,scale= rowsum.weightsMatrix))
datXscaledweighted= as.matrix(datXscaled* weightsMatrix)
# this corresponds to formula 25 in Klemera et al 2006
y.true2=as.numeric(apply(datXscaledweighted,1,sum,na.rm=TRUE))
 } #end of if (skipMissingVariables )


# the following is different from Klemera in that it has an absolute value
r.characteristic=sum(rVector^2/sqrt(1-rVector^2) )/sum(abs(rVector)/sqrt(1-rVector^2) )
no.missing=sum(apply(  as.matrix( is.na(datX)),1,sum))
if (sum(no.missing)>0) {warning("The input datX contains missing values.\n
I recommend you impute missing values in datX before running this function. ") 
} # end of if (sum(no.missing)>0)
# formula 37 from Klemera
SsquaredBE=var( y.true2-y,na.rm=TRUE) -(1- r.characteristic^2)/r.characteristic^2*var(y,na.rm=TRUE)/no.variables
# this corresponds to formula 34 in Klemera
y.true3=(as.numeric( as.matrix(datXscaled)%*% weights0)+y/SsquaredBE )/( sum(weights0)+ 1/SsquaredBE)
y.true3[is.na(y.true3) ]=y[is.na(y.true3)]
} # end of if (no.observations>1 ) 
SD.ytrue2=sqrt(1-r.characteristic^2)/r.characteristic*sqrt(var(y,na.rm=TRUE)/no.variables)
# now formula 42
SD.ytrue3=SD.ytrue2/sqrt(1+SD.ytrue2^2/SsquaredBE )
} # end of if (sum(restNonMissingY,na.rm=TRUE) >3 )
datEstimates=data.frame(y, y.true1,y.true2,y.true3)

if (!is.null(datXtest)){
datXtest=as.matrix(datXtest)
no.variablestest=dim(as.matrix(datXtest))[[2]]
if (no.variablestest != no.variables) {stop("the number of variables in the test data is not the same as in the training data")}

y.true1test=rep(NA, length(y) )
y.true2test=rep(NA, length(y) )
y.true3test=rep(NA, length(y) )
restNonMissingY= !is.na(y)
if (sum(restNonMissingY,na.rm=TRUE) >3 ){
datXtestscaled=scale(datXtest,center=datCoef[,1],scale=datCoef[,2] )
y.true1test=as.numeric(apply(as.matrix(datXtestscaled),1,mean) )
y.true2test=as.numeric( as.matrix(datXtestscaled)%*%weights)
if (skipMissingVariables ) {
y.true1test= as.numeric(apply(as.matrix(datXtestscaled),1,mean,na.rm=TRUE))
weightsMatrixtest=matrix(weights,byrow=TRUE,nrow=dim(as.matrix(datXtestscaled))[[1]],ncol=length(weights) )
weightsMatrixtest[is.na(datXtestscaled)]=0
rowsum.weightsMatrixtest=apply(as.matrix(weightsMatrixtest),1,sum)
weightsMatrixtest=t(scale(t(as.matrix(weightsMatrixtest)),center=F,scale= rowsum.weightsMatrixtest))
datXscaledweightedtest= as.matrix(datXtestscaled* weightsMatrixtest)
# this corresponds to formula 25 in Klemera et al 2006
y.true2test=as.numeric(apply(datXscaledweightedtest,1,sum,na.rm=TRUE))
 } #end of if (skipMissingVariables )
} # end of if (sum(restNonMissingY,na.rm=TRUE) >3 )


datEstimatestest=data.frame(y.true1= y.true1test,y.true2= y.true2test)

} # end of if (!is.null(datXtest))

if ( LeaveOneOut.CV  ) {
y.true1test.LOO=rep(NA,no.observations)
y.true2test.LOO=rep(NA,no.observations)
y.lmLOO= rep(NA,no.observations)
for ( i in 1:no.observations ){
rm(datCoef); rm(corX);

datX.LOO=datX[-i,]
datXtest.LOO= matrix(datX[i,],nrow=1)
y.LOO=y[-i]
no.variables=dim(as.matrix(datX.LOO))[[2]]
no.observations=dim(as.matrix(datX.LOO))[[1]]

if (no.observations==1 ) {
warning("When dealing with leave one out cross validation, there is only 1 observations in the training data")}

if (no.observations>1 ) {

if (addLinearModel) {
lmLOO=lm(y.LOO~., data=data.frame(datX.LOO),na.action=na.exclude)
y.lmLOO[i]= sum(datXtest.LOO*lmLOO$coeff[-1])+lmLOO$coeff[[1]]
}

corX = parse(text = paste(corFnc, "(datX.LOO,y.LOO ", 
                prepComma(corOptions), ")"))
rVector= as.numeric(eval(corX))
datCoef=t(coef(lm(datX.LOO~ y.LOO,na.action="na.exclude")))
datX.LOOscaled=scale(datX.LOO,center=datCoef[,1],scale=datCoef[,2] )
weights0=rVector^2/(1-rVector^2)
weights=weights0/sum(weights0)
datXtest.LOOscaled=(datXtest.LOO-datCoef[,1])/datCoef[,2] 
y.true1test.LOO[i]= mean(datXtest.LOOscaled)
y.true2test.LOO[i]=sum(datXtest.LOOscaled*weights)
if (skipMissingVariables ) {
y.true1test.LOO[i]= mean(datXtest.LOOscaled,na.rm=TRUE)
weightsMatrixLOO= weights
weightsMatrixLOO[is.na(datXtest.LOOscaled)]=0
rowsum.weightsMatrixLOO=sum(weightsMatrixLOO)
weightsMatrixLOO=weightsMatrixLOO/rowsum.weightsMatrixLOO
datXscaledweightedLOO= datXtest.LOOscaled* weightsMatrixLOO
y.true2test.LOO[i]=sum(datXscaledweightedLOO,na.rm=TRUE)
 } #end of if (skipMissingVariables )
} # end of for loop
} # end of if (no.observations>1 ) 
datEstimates.LeaveOneOut.CV=data.frame(y.true1= y.true1test.LOO,
y.true2= y.true2test.LOO)
} # end of if ( LeaveOneOut.CV  ) 

if (addLinearModel) {
y.lmTest=rep(NA, dim(as.matrix(datXtest))[[1]] )
restNonMissingY= !is.na(y)
if (sum(restNonMissingY,na.rm=TRUE) >3 ){
lm1=lm(y~., data=data.frame(datX),na.action=na.exclude)
y.lmTraining=predict(lm1)
y.lmTraining=predict(lm1)
if( !is.null(datXtest)) {y.lmTest=predict(lm1,newdata=data.frame(datXtest))}
}
}
if ( !is.null(datXtest) & LeaveOneOut.CV & !addLinearModel ) {out=  list( datEstimates=datEstimates, datEstimatestest= datEstimatestest,
datEstimates.LeaveOneOut.CV= datEstimates.LeaveOneOut.CV, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( !is.null(datXtest) & !LeaveOneOut.CV & !addLinearModel) {out=  list( datEstimates=datEstimates, datEstimatestest= datEstimatestest, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( is.null(datXtest) & LeaveOneOut.CV & !addLinearModel) {out=  list( datEstimates=datEstimates, datEstimates.LeaveOneOut.CV= datEstimates.LeaveOneOut.CV, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( is.null(datXtest) & !LeaveOneOut.CV & !addLinearModel) {out=  list( datEstimates=datEstimates, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }

if ( !is.null(datXtest) & LeaveOneOut.CV & addLinearModel ) {out=  list( datEstimates=data.frame(datEstimates, y.lm=y.lmTraining), datEstimatestest= data.frame(datEstimatestest, y.lm=y.lmTest),
datEstimates.LeaveOneOut.CV= data.frame(datEstimates.LeaveOneOut.CV,y.lm=y.lmLOO), SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }

if ( !is.null(datXtest) & !LeaveOneOut.CV & addLinearModel) {out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining), datEstimatestest= data.frame(datEstimatestest,y.lm=y.lmTest), SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }

if ( is.null(datXtest) & LeaveOneOut.CV & addLinearModel) {out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining), datEstimates.LeaveOneOut.CV= data.frame(datEstimates.LeaveOneOut.CV,y.lm=y.lmLOO),
 SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( is.null(datXtest) & !LeaveOneOut.CV & addLinearModel) {out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining), SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
out
} # end of function


