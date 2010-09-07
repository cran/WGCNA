.absMax <- function(datIn){
	datIn = abs(datIn)
	keep = which.max(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.absMin <- function(datIn){
	datIn = abs(datIn)
	keep = which.min(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.Max <- function(datIn){
# Note that this is a matrix version of the "max" function
	keep = which.max(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.Min <- function(datIn){
# Note that this is a matrix version of the "min" function
	keep = which.min(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.maxRowVariance <- function(datIn){
	sds  = apply(datIn,1,function(x) return(var(x,na.rm=TRUE)))
	keep = which.max(sds)[1]
	return(as.numeric(datIn[keep,]))
}

.selectFewestMissing <- function(datET, probeVector, symbolVector, omitGenes, omitPercent=90){
## For each gene, select the gene with the fewest missing probes, and return the results.
#   If there is a tie, keep all probes involved in the tie.
#   The main part of this funciton is run only if omitGenes=TRUE
	
	# First, return datET if there is no missing data, otherwise run the function
	if (sum(is.na(datET))==0) return(datET)
	
	# Set up the variables.
	names(symbolVector) = probeVector
	probes              = rownames(datET)
	genes               = symbolVector[probes]
	keepGenes           = rep(TRUE,length(probes))
	tGenes              = table(genes)
	checkGenes          = sort(names(tGenes)[tGenes>1])
	missingData         = rowSums(is.na(datET))
	
	# Omit all probes with at least omitPercent genes missing
	datET = datET[missingData<(omitPercent*dim(datET)[2]/100),]
	
	# Omit relevant genes and return results
	if (omitGenes)
		for (g in checkGenes){
			gn            = (genes==g)
			keepGenes[gn] = (missingData[gn] == min(missingData[gn]))
		}
	return(datET[keepGenes,])
}



# ----------------- Main Function ------------------- #

collapseRows <- function(datET, symbolVector, probeVector, method="maxRowVariance", connectivityBasedCollapsing=TRUE,	
	methodFunction=NULL, connectivityPower=1, selectFewestMissing=TRUE){

# datET is an expression matrix with rows=genes and cols=samples

## Test to make sure the variables are the right length.
#     if not, fix it if possible, or return 0 if not possible
	probeVector  = as.character(probeVector)
	symbolVector = as.character(symbolVector)
	rnDat = rownames(datET)
	if (length(probeVector)!=length(symbolVector)){
		write("Error: symbolVector and probeVector not the same length... exiting.","")
		return(0)
	}
	if ((is.null(rnDat))&(dim(datET)[1]==length(probeVector))){
		write("Warning: *datET* does not have row names.  Assigning *probeVector* as row names.","")
		rnDat <- rownames(datET) <- probeVector
	}
	if (is.null(rnDat)){
		write("Error: *datET* does not have row names and length of *probeVector*...","")
		write("... is not the same as # rows in *datET*... exiting.","")
		return(0)
	}
	if (sum(is.element(rnDat,probeVector))!=length(probeVector)){
		write("Warning: row names of input data and probes not identical...","")
		write("... Attempting to proceed anyway. Check results carefully.","")
	}
		
## For each gene, select the gene with the fewest missing probes (if selectFewestMissing==TRUE)
##  Also, remove all probes with more than 90% missing data
	
	datET = .selectFewestMissing(datET, probeVector, symbolVector, selectFewestMissing)
	rnDat = rownames(datET)
		
##   If method="function", use the function "methodFunction" as a way of combining genes
#    Alternatively, use one of the built-in functions 
#    Note: methodFunction must be a function that takes a vector of numbers as input and
#     outputs a single number. This function will return(0) or crash otherwise.

        imethod = match(method, c("function","Max","maxRowVariance","Min","absMin","absMax"));
        
	if (is.na(imethod)) {
		write("Error: entered method is not a legal option.  Please enter *maxRowVariance*,","")
		write(" *Max*, *Min*, *absMax*, *absMin*, or *function* for a user-defined function.","")
		return(0)
	}
        if (imethod > 1) method = spaste(".", method);
	if(method=="function") method = methodFunction
	if((!is.function(methodFunction))&(!is.null(methodFunction))){
		write("Error: *methodFunction* must be a function... please read the help file","")
		return(0)
	}
	method = match.fun(method)
		
## Format the variables for use by this function
	probeVector[is.na(probeVector)] = symbolVector[is.na(probeVector)] # Use gene if probe is missing
	rownames(datET)[is.na(rnDat)]   = symbolVector[is.na(rnDat)]
	remove       = (is.na(probeVector))|(is.na(symbolVector)) # Omit if both gene and probe are missing
	probeVector  = probeVector[!remove];
	symbolVector = symbolVector[!remove];
	names(symbolVector) = probeVector
	probeVector = sort(intersect(rnDat,probeVector))
	if (length(probeVector)==0){
		write("Error: none of the *datET* rownames are in *probeVector*...","")
		write("... please add rownames and try again... exiting.","")
		return(0)
	}
	symbolVector = symbolVector[probeVector]
	datET  = as.matrix(datET)
	datET  = datET[probeVector,]
	probes = rownames(datET)
	genes  = symbolVector[probes]
	tGenes = table(genes)
	datETOut=matrix(0,nrow=length(tGenes),ncol=length(colnames(datET)))
	colnames(datETOut) = colnames(datET)
	rownames(datETOut) = sort(names(tGenes))
	probesOut = rownames(datETOut)
	names(probesOut) = probesOut
	
##  If !is.null(connectivityPower), default to the connectivity method with power=method
#      Collapse genes with multiple probe sets together using the following algorthim:
#      1) If there is one ps/g = keep
#      2) If there are 2 ps/g = (use "method" or "methodFunction")
#      3) If there are 3+ ps/g = take the max connectivity
#   Otherwise, use "method" if there are 3+ ps/g as well. 
	if(!is.null(connectivityPower)){
		if(!is.numeric(connectivityPower)){
			write("Error: if entered, connectivityPower must be numeric... exiting.","")
			return(0)
		}
		if(connectivityPower<=0){
			write("Warning: connectivityPower must be >= 0.  Defaulting to a power of 2.","")
			connectivityPower=2
		}
		if(dim(datET)[2]<=5){
			write("Warning: 5 or fewer samples, this method of probe collapse is unreliable...","")
			write("...Running anyway, but we suggest trying another method (for example, *mean*).","")
		}
	}
	
# Actually run the collapse now!!!
	if (!is.null(methodFunction))
		write("Comment: make sure methodFunction takes a matrix as input.","")
	ones = sort(names(tGenes)[tGenes==1])
	if(connectivityBasedCollapsing){
		twos = sort(names(tGenes)[tGenes==2]) # use "method" and connectivity
		more = sort(names(tGenes)[tGenes>2])
	} else { 
		twos = sort(names(tGenes)[tGenes>1]) # only use "method"
		more = character(0)
	}
	for (g in ones){
		datETOut[g,] = as.numeric(datET[probes[genes==g],])
		probesOut[g] = probes[genes==g]
	}
	for (g in twos){
		datETTmp = datET[probes[genes==g],]
		datETOut[g,] = as.numeric(method(datETTmp))
		whichTest    = apply(datETTmp,1,function(x) return(sum(x==datETOut[g,])))
		probesOut[g] = (names(whichTest)[whichTest==max(whichTest)])[1]
	}
	for (g in more){
		datETTmp = datET[probes[genes==g],]
		adj = (0.5+0.5*cor(t(datETTmp),use="p"))^connectivityPower
		datETOut[g,] = as.numeric(datETTmp[which.max(rowSums(adj,na.rm=TRUE)),])
		whichTest    = apply(datETTmp,1,function(x) return(sum(x==datETOut[g,])))
		probesOut[g] = (names(whichTest)[whichTest==max(whichTest)])[1]
	}
	if (!is.null(methodFunction))
		write("...Ignore previous comment.  Function completed properly!","")
		
# Retreive the information about which probes were saved, and include that information
#   as part of the output
	out2 = cbind(rownames(datETOut),probesOut)
	colnames(out2) = c("GeneSymbol","SelectedProbes")
	out3 = is.element(rownames(datET),probesOut)
	names(out3) = rownames(datET)
	output = list(expression = datETOut, symbol2probe = out2, originalProbes = out3)
	return(output)
		
# End of function
} 
