collapseRows <- function(datET, symbolVector, probeVector, method=max, connectivityPower=NULL){

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
	
##   Use the function "method" as a way of combining genes
#    Note: method must be a function that takes a vector of numbers as input and
#     outputs a single number. This function will return(0) or crash otherwise.
	if(!is.function(method)){
		write("Error: *method* must be function (default=mean)... please read help file","")
		return(0)
	}
	method = match.fun(method)
	
# This internal function removes the NAs before passing the variables onto 
#   the "method" function.  If there are only NAs, NA is returned.
	collapseMethod = function(x){
		xOut = x[!is.na(x)]
		if(length(xOut)==0) return(NA)
		return(method(xOut))
	} # End internal function
	
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
	datET = as.matrix(datET)
	datET = datET[probeVector,]
	probes = rownames(datET)
	genes  = symbolVector[probes]
	tGenes = table(genes)
	datETOut=matrix(0,nrow=length(tGenes),ncol=length(colnames(datET)))
	colnames(datETOut) = colnames(datET)
	rownames(datETOut) = sort(names(tGenes))
	
##  If !is.null(connectivityPower), default to the connectivity method with power=method
#      Collapse genes with multiple probe sets together using the following algorthim:
#      1) If there is one ps/g = keep
#      2) If there are 2 ps/g = (use "method")
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
	write("Comment: make sure *method* function takes a single numeric vector as input.","")
	ones = sort(names(tGenes)[tGenes==1])
	if(!is.null(connectivityPower)){
		twos = sort(names(tGenes)[tGenes==2]) # use "method" and connectivity
		more = sort(names(tGenes)[tGenes>2])
	} else { 
		twos = sort(names(tGenes)[tGenes>1]) # only use "method"
		more = character(0)
	}
	for (g in ones)
		datETOut[g,] = as.numeric(datET[probes[genes==g],])
	for (g in twos){
		datETTmp = datET[probes[genes==g],]
		datETOut[g,] = as.numeric(apply(datETTmp,2,collapseMethod))
	}
	for (g in more){
		datETTmp = datET[probes[genes==g],]
		adj = cor(t(datETTmp))^connectivityPower
		datETOut[g,] = as.numeric(datETTmp[which.max(rowSums(adj,na.rm=TRUE)),])
	}
	write("...Ignore previous comment.  Function completed properly!","")
	return(datETOut)
# End of function
} 
