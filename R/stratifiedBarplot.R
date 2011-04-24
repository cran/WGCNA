# This function barplots data across two splitting parameters
stratifiedBarplot = function (expAll, groups, split, subset, 
genes=NA, scale="N", graph=TRUE, las1=2, cex1=1.5, ...){

## Code to take care of array formatting
	expAll = t(expAll)
	if (length(subset)>1){
		subsetName = subset[1]; 
		subset = subset[2:length(subset)]
	} else { subsetName=subset }
	groupNames = as.character(names(.tableOrd(groups)))
	splitNames = as.character(names(.tableOrd(split))) 
	if (is.na(genes)[1]) genes = rownames(expAll)
	keep   = !is.na(groups)
	groups = groups[keep]
	split  = split[keep]
	expAll = expAll[,keep]
	scale  = substr(scale,1,1)
	
## Collect and scale the expression data
	expSubset = expAll[is.element(genes,subset),]
	if(length(subset)>1){
		if(scale=="A")  expSubset = t(apply(expSubset,1,function(x) return(x/mean(x))))
		if(scale=="Z")  expSubset = t(apply(expSubset,1,function(x) return((x-mean(x))/sd(x))))
		if(scale=="H")  {
			AdjMat = adjacency(t(expSubset),type="signed",power=2)
			diag(AdjMat) = 0
			Degree = rowSums(AdjMat)
			keep   = which(Degree == max(Degree))
			expSubset = expSubset[keep,]	
		}
		if(scale=="M")  {
			me = moduleEigengenes(as.matrix(t(expSubset)), rep("blue",dim(expSubset)[1]))
			expSubset = me$eigengenes$MEblue
		} 
		expSubset = rbind(expSubset,expSubset)
		expSubset = apply(expSubset,2,mean)
	}
	
## Now average the data and output it in a meaningful way
	exp <- std <- matrix(0,ncol=length(splitNames),nrow=length(groupNames))
	splitPvals = rep(1,length(splitNames))
	names(splitPvals) = splitNames
	groupPvals = rep(1,length(groupNames))
	names(groupPvals) = groupNames
	for (c in 1:length(splitNames)){
		expTmp = expSubset[split==splitNames[c]]
		grpTmp = groups[split==splitNames[c]]
		splitPvals[c] = kruskal.test(expTmp,as.factor(grpTmp))$p.value
		for (r in 1:length(groupNames)){
			exp[r,c] = mean(expSubset[(groups==groupNames[r])&(split==splitNames[c])])
			std[r,c] = sd(expSubset[(groups==groupNames[r])&(split==splitNames[c])])
			if(c==1){
				expTmp = expSubset[groups==groupNames[r]]
				splTmp = split[groups==groupNames[r]]
				groupPvals[r] = kruskal.test(expTmp,as.factor(splTmp))$p.value	
			}
		}
	}
	colnames(exp) <- colnames(std) <- splitNames
	rownames(exp) <- rownames(std) <- groupNames
	
## Now plot the results, if requested
	if(graph){
		ylim = c(min(0,min(min(exp-std))),max(max(exp+std))*(1+0.08*length(groupNames)))
		barplot(exp, beside=TRUE, legend.text=TRUE, main=subsetName, las=las1,
				ylim=ylim, cex.axis=cex1, cex.names=cex1,  ...)
		.err.bp(exp,std,TRUE)
	}
	
# Now collect the output and return it.
	out = list(splitGroupMeans = exp, splitGroupSDs = std, 
			   splitPvals = splitPvals, groupPvals = groupPvals)
	return(out)
}
# --------------------------------


.err.bp<-function(daten,error,two.side=F){
## This function was written by Steve Horvath
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)
	if(!is.numeric(daten)) {
	stop("All arguments must be numeric")}
	if(is.vector(daten)){ 
		xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
	}else{
		if (is.matrix(daten)){
			xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
							   dim=c(1,length(daten))))+0:(length(daten)-1)+.5
		}else{
		stop("First argument must either be a vector or a matrix") }
	}
	MW<-0.25*(max(xval)/length(xval)) 
	ERR1<-daten+error 
	ERR2<-daten-error
	for(i in 1:length(daten)){
		segments(xval[i],daten[i],xval[i],ERR1[i])
		segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
		if(two.side){
			segments(xval[i],daten[i],xval[i],ERR2[i])
			segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
		} 
	} 
} 


.tableOrd = function (input){
## This is the same as the "table" function but retains the order
	
## This internal function collects the order
	tableOrd2 = function(input, output=NULL){
		input  = input[!is.na(input)]
		if (length(input)==0) return (output)
		outTmp = input[1] 
		output = c(output, outTmp)
		input  = input[input!=outTmp]
		output = tableOrd2(input, output)
		return(output)
	}
	
## Get the results
	tableOut   = table(input)
	tableOrder = tableOrd2(input)
	return(tableOut[tableOrder])
}  
