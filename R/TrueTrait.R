TrueTrait=function(datX, y,datXtest=NULL, corFnc = "cor", corOptions = "use = 'pairwise.complete.obs'", LeaveOneOut.CV=FALSE,skipMissingVariables=TRUE,addLinearModel=FALSE, Strata=NULL){
datX=as.matrix(datX)
no.variables=dim(as.matrix(datX))[[2]]
r.characteristic=NA
datVariableInfo=data.frame(matrix(NA, nrow=no.variables, ncol=4))
names(datVariableInfo)=c("Variable","center", "scale", "weights.y.true2")
if ( is.null(colnames(datX))){datVariableInfo$Variable=1:no.variables} else { datVariableInfo$Variable=colnames(datX)}
no.observations=dim(as.matrix(datX))[[1]]
if (no.observations !=  length(y) ) {stop("The number of rows of datX does not correspond to the length of y. Consider transposing your input matrix or use a different input.") }
if (no.observations==1 ) {
warning("Only 1 observations, i.e. the length of y is 1. The function cannot be used. For your convenience, the estimated true values will be set to the input value of y.")
y.true1=y; y.true2=y;  y.true3=y}
if (no.observations>1 ) {
vary=var(y,na.rm=TRUE)
meany=mean(y,na.rm=TRUE)
if (vary==0 | is.na(vary) ) { 
y.true1=rep(meany ,length(y))  ; y.true2= rep(meany ,length(y)) ;  y.true3= rep(meany ,length(y))     }
if (vary> 0 & !is.na(vary) ) { 
corX = parse(text = paste(corFnc, "(datX,y ",prepComma(corOptions), ")"))
rVector= as.numeric(eval(corX))
datCoef=t(coef(lm(datX~ y,na.action="na.exclude")))
datVariableInfo$center=datCoef[,1]
datVariableInfo$scale=datCoef[,2]
datXscaled=scale(datX,center=datCoef[,1],scale=datCoef[,2] )
weights0=rVector^2/(1-rVector^2)
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
SsquaredBE=var( y.true2-y) -(1- r.characteristic^2)/r.characteristic^2*var(y)/no.variables
# this corresponds to formula 34 in Klemera
y.true3=(as.numeric( as.matrix(datXscaled)%*% weights0)+y/SsquaredBE )/( sum(weights0)+ 1/SsquaredBE)
y.true3[is.na(y.true3) ]=y[is.na(y.true3)]
} # end of if (vary> 0 & !is.na(vary) ) 
} # end of if (no.observations>1 ) 

SD.ytrue2=NA
SD.ytrue3=NA
if (vary> 0 & !is.na(vary) & r.characteristic>0 & !is.na(r.characteristic) ) { 
SD.ytrue2=sqrt(1-r.characteristic^2)/r.characteristic*sqrt(var(y)/no.variables)
# now formula 42
SD.ytrue3=SD.ytrue2/sqrt(1+SD.ytrue2^2/SsquaredBE )
} # end of if if (vary> 0 & !is.na(vary)

datEstimates=data.frame(y.true1,y.true2)

if (!is.null(datXtest)){
datXtest=as.matrix(datXtest)
no.variablestest=dim(as.matrix(datXtest))[[2]]
if (no.variablestest != no.variables) {stop("the number of variables in the test data is not the same as in the training data")}

if (vary==0 | is.na(vary) ) { 
no.testset= dim( as.matrix(datXtest) )[[1]]
y.true1test=rep(meany, no.testset)  ; y.true2test= rep(meany, no.testset)   ;  y.true3test= rep(meany, no.testset)       }
if (vary> 0 & !is.na(vary) ) { 
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
} # end of if (vary> 0 & !is.na(vary) ) 

datEstimatestest=data.frame(y.true1= y.true1test,y.true2= y.true2test)

} # end of if (!is.null(datXtest))

if ( LeaveOneOut.CV  ) {
y.true1test.LOO=rep(NA,no.observations)
y.true2test.LOO=rep(NA,no.observations)
y.lmLOO= rep(NA,no.observations)
if (vary==0 | is.na(vary) ) { 
y.true1test.LOO=rep(meany, no.observations)  ; y.true2test.LOO= rep(meany, no.observations)   ;  y.true3test.LOO= rep(meany, no.observations)     
y.lmLOO= rep(meany,no.observations) }
if (vary> 0 & !is.na(vary) ) { 
for ( i in 1:no.observations ){
 rm(datCoef); rm(corX)
datX.LOO=datX[-i,]
datXtest.LOO= matrix(datX[i,],nrow=1)
y.LOO=y[-i]
vary.LOO=var(y.LOO,na.rm=TRUE)
meany.LOO= mean(y.LOO,na.rm=TRUE)
no.variables=dim(as.matrix(datX.LOO))[[2]]
no.observations=dim(as.matrix(datX.LOO))[[1]]

if (no.observations==1 ) {
warning("When dealing with leave one out cross validation, there is only 1 observations in the training data")}

if (no.observations>1 ) {

if (addLinearModel) {
if ( vary.LOO==0 | is.na(vary.LOO) ) {y.lmLOO[i]=meany.LOO} else {
lmLOO=lm(y.LOO~., data=data.frame(datX.LOO),na.action=na.exclude)
y.lmLOO[i]= sum(datXtest.LOO*lmLOO$coeff[-1])+lmLOO$coeff[[1]]
}
}

corX = parse(text = paste(corFnc, "(datX.LOO,y.LOO ", 
                WGCNA:::prepComma(corOptions), ")"))
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
} # end of if (vary> 0 & !is.na(vary) ) 
} # end of if (no.observations>1 ) 

datEstimates.LeaveOneOut.CV=data.frame(y.true1= y.true1test.LOO,
y.true2= y.true2test.LOO)
} # end of if ( LeaveOneOut.CV  ) 

if (addLinearModel) {

if (vary==0 | is.na(vary) ) { 
y.lmTraining=rep(meany,length(y) )
if( !is.null(datXtest)) {y.lmTest=rep(meany, dim(as.matrix(datXtest))[[1]] )  }
}
if (vary> 0 & !is.na(vary) ) { 
lm1=lm(y~., data=data.frame(datX),na.action=na.exclude)
y.lmTraining=predict(lm1)
if( !is.null(datXtest)) {y.lmTest=predict(lm1,newdata=data.frame(datXtest))}
} # end of if (vary> 0 & !is.na(vary) ) 
} # if (addLinearModel) 



if ( !is.null(datXtest) & LeaveOneOut.CV & !addLinearModel ) {out=  list( datEstimates=data.frame(datEstimates, ,y.true3), datEstimatestest= datEstimatestest,
datEstimates.LeaveOneOut.CV= datEstimates.LeaveOneOut.CV, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( !is.null(datXtest) & !LeaveOneOut.CV & !addLinearModel) {out=  list( datEstimates=datEstimates, datEstimatestest= datEstimatestest, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( is.null(datXtest) & LeaveOneOut.CV & !addLinearModel) {out=  list( datEstimates=datEstimates, datEstimates.LeaveOneOut.CV= datEstimates.LeaveOneOut.CV, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( is.null(datXtest) & !LeaveOneOut.CV & !addLinearModel) {out=  list( datEstimates=datEstimates, SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }

if ( !is.null(datXtest) & LeaveOneOut.CV & addLinearModel ) {out=  list( datEstimates=data.frame(datEstimates, y.lm=y.lmTraining, y.true3), datEstimatestest= data.frame(datEstimatestest, y.lm=y.lmTest),
datEstimates.LeaveOneOut.CV= data.frame(datEstimates.LeaveOneOut.CV,y.lm=y.lmLOO), SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }


if ( !is.null(datXtest) & !LeaveOneOut.CV & addLinearModel) {out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining, y.true3), datEstimatestest= data.frame(datEstimatestest,y.lm=y.lmTest), SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }

if ( is.null(datXtest) & LeaveOneOut.CV & addLinearModel) {out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining, y.true3), datEstimates.LeaveOneOut.CV= data.frame(datEstimates.LeaveOneOut.CV,y.lm=y.lmLOO),
 SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
if ( is.null(datXtest) & !LeaveOneOut.CV & addLinearModel) {out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining,y.true3), SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }



if (   !is.null(Strata)  )  {
if ( dim(datX)[[1]]  != length(Strata) ) stop("The length of the Strata vector does not equal the number of observations (rows of datX)")


uniqueStrata=unique(Strata)
outByStratum=out
if ( !is.null(datXtest) ) outByStratum$datEstimatestest=NULL
names(outByStratum)=paste(names(outByStratum),"ByStratum",sep="")
if (  length(uniqueStrata) ==1 ) { outByStratum=out; warning("There is only one stratum. No need to specify the variable Strata") } else {
outByStratum$SD.ytrue2ByStratum=rep(NA, length(uniqueStrata) )
outByStratum$SD.ytrue3ByStratum= rep(NA, length(uniqueStrata) )
outByStratum$datVariableInfoByStratum=list(NULL)


for (i in 1:length(uniqueStrata) ) {
whichStratum=uniqueStrata[i] 
#printFlush(i)
#printFlush(whichStratum)
if ( is.na(whichStratum) ) {rest2=is.na(Strata) } else { 
rest2=Strata==whichStratum & !is.na(Strata)  }
if (whichStratum=="exclude" | is.na(whichStratum) ) {
if (is.element("datEstimates",  names(out)) ) {outByStratum$datEstimatesByStratum[ rest2,]=
matrix(NA,nrow=nrow(outByStratum$datEstimatesByStratum[ rest2,]) ,ncol=ncol(outByStratum$datEstimatesByStratum[ rest2,]) )}
if (is.element("datEstimates.LeaveOneOut.CV",  names(out)) ) {outByStratum$datEstimates.LeaveOneOut.CVByStratum[ rest2,]=
matrix(NA,nrow=nrow(outByStratum$datEstimates.LeaveOneOut.CV[ rest2,]) ,ncol=ncol(out$datEstimates.LeaveOneOut.CV[ rest2,]) )}
outByStratum$SD.ytrue2ByStratum[i]=NA
outByStratum$SD.ytrue3ByStratum[i]=NA
outByStratum$datVariableInfoByStratum[[i]]=matrix(NA,nrow=nrow(outByStratum$datVariableInfoByStratum[[1]]), ncol=ncol(outByStratum$datVariableInfoByStratum[[1]])) 
} # end of  if (whichStratum=="exclude" | is.na(whichStratum) )
if (whichStratum !="exclude" & !is.na(whichStratum) ) {
y2=y[rest2]
datX2= datX[rest2,]
outStratum=TrueTrait(datX=datX2, y=y2, datXtest=NULL , LeaveOneOut.CV=LeaveOneOut.CV, addLinearModel=addLinearModel,corFnc = corFnc, corOptions = corOptions, skipMissingVariables=skipMissingVariables,Strata=NULL)
outByStratum$datEstimatesByStratum[ rest2,]=outStratum$datEstimates
if (is.element("datEstimates.LeaveOneOut.CV",  names(out)) ) {outByStratum$datEstimates.LeaveOneOut.CVByStratum[ rest2,]=outStratum$datEstimates.LeaveOneOut.CV}
outByStratum$SD.ytrue2ByStratum[i]=outStratum$SD.ytrue2
outByStratum$SD.ytrue3ByStratum[i]=outStratum$SD.ytrue3
outByStratum$datVariableInfoByStratum[[i]]=outStratum$datVariableInfo
} # end of if (whichStratum !="exclude" )
} # end of for (i in 1:length(uniqueStrata) )
}  # end of if (length(uniqueStrata)>1) 
out=append( out, outByStratum)
} # end of if (   !is.null(Strata)  )  
out
} # end of function





