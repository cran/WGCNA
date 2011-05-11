.filterSimilarPS <- function(datOut, rowGroup, rowID, thresh=0.8){
## Collapse groups (ie. genes) with multiple ids (ie. probes) together using the following algorthim:
# 1) If there is one id/group = keep
# 2) If there are 2 ids/group = take the maximum mean expression, if their correlation is > thresh
# 3) If there are 3+ ids/group = iteratively repeat (2) for the id with the highest 
#	correlation until all ids remaining have correlation < thresh for each group
# datOut is an expression matrix with rows=ids (NOT group) and cols=samples
# rowGroup and rowID are vectors of corresponding group and id names for the rows in datOut
#	(note: all ids in datOut only need to be a subset of these vectors, not necessarily identical)
# thresh is the Pearson correlation threshold to combine probes of similar expression.
	
	names(rowGroup) = rowID
	ids    = rownames(datOut);  
	idsIn  = ids;   # For later
	group  = rowGroup[ids]
	tGroup = table(group)
	twos   = sort(names(tGroup)[tGroup==2])
	more   = sort(names(tGroup)[tGroup>2])
	len    = dim(datOut)[2]
	
	testTwoAndCombine <- function (datIn, datTmp, thresh){
	# Internal function for removing one of two genes if they have high enough correlation
		if (cor(as.numeric(datTmp[1,]),as.numeric(datTmp[2,]))>thresh){
			rMean = rowMeans(datTmp)
			omit  = as.numeric(which(rMean==min(rMean)))
			datIn = datIn[rownames(datIn)!=rownames(datTmp)[omit],]
		} 
		return(datIn)
	} # End internal function "testTwoAndCombine"
	
	for (g in twos){
		datTmp = datOut[ids[group==g],]
		datOut = testTwoAndCombine(datOut, datTmp, thresh)
	}; write("Done combining genes with 2 probes!","")
	for (g in more){
		go = TRUE
		while(go){
			ids = rownames(datOut)
			group  = rowGroup[ids]
			datTmp = datOut[ids[group==g],]
			corDat = cor(t(datTmp)); 
			if(length(datTmp)==len) { go=FALSE
			} else {
				diag(corDat)=-2
				if (max(corDat)<thresh) { go=FALSE
				} else {
					w  = which(corDat==max(corDat))[1]
					l  = dim(datTmp)[1]
					wI = c(((w-1)%%l)+1,ceiling(w/l))
					datTmp = datTmp[wI,]
					datOut = testTwoAndCombine(datOut, datTmp, thresh)
				}
			}
		}
	}; write("Done combining genes with 3+ probes!","")
	idsNew      = rownames(datOut)
	groupNew    = rowGroup[idsNew]
	out2        = cbind(groupNew,idsNew)
	selectedRow = is.element(idsIn,idsNew)
	names(selectedRow) = ids
	output      = list(datETcollapsed = datOut, group2row = out2, selectedRow = selectedRow)
	return(output)
}

.absMaxMean <- function(datIn){
	datIn = abs(datIn)
	keep = which.max(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.absMinMean <- function(datIn){
	datIn = abs(datIn)
	keep = which.min(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.MaxMean <- function(datIn){
# Note that this is a matrix version of the "max" function
	keep = which.max(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.MinMean <- function(datIn){
# Note that this is a matrix version of the "min" function
	keep = which.min(rowSums(datIn,na.rm=TRUE))[1]
	return(as.numeric(datIn[keep,]))
}

.maxRowVariance <- function(datIn){
	sds  = apply(datIn,1,function(x) return(var(x,na.rm=TRUE)))
	keep = which.max(sds)[1]
	return(as.numeric(datIn[keep,]))
}

.Average <- function(datIn){
	return(as.numeric(colMeans(datIn)))
}

.selectFewestMissing <- function(datET, rowID, rowGroup, omitGroups, omitPercent=90){
## For each gene, select the gene with the fewest missing probes, and return the results.
#   If there is a tie, keep all probes involved in the tie.
#   The main part of this function is run only if omitGroups=TRUE
	
	# First, return datET if there is no missing data, otherwise run the function
	if (sum(is.na(datET))==0) return(rep(TRUE,nrow(datET)));
	
	# Set up the variables.
	names(rowGroup)     = rowID
	probes              = dimnames(datET)[[1]]
	genes               = rowGroup[probes]
	keepGenes           = rep(TRUE,length(probes))
	tGenes              = table(genes)
	checkGenes          = sort(names(tGenes)[tGenes>1])
	missingData         = rowSums(is.na(datET))
	
	# Omit all probes with at least omitPercent genes missing
        keep = missingData<(omitPercent*dim(datET)[2]/100);
	
	# Omit relevant genes and return results
	if (omitGroups)
		for (g in checkGenes){
			gn            = (genes==g)
			keepGenes[gn] = (missingData[gn] == min(missingData[gn]))
		}

        keep = keep & keepGenes;
	return (keep);
}

# ----------------- Main Function ------------------- #

collapseRows <- function(datET, rowGroup, rowID, method="MaxMean", connectivityBasedCollapsing=FALSE,	
	methodFunction=NULL, connectivityPower=1, selectFewestMissing=TRUE, thresholdCombine=NA)
{

	# datET = as.matrix(as.data.frame(datET));
	methodAverage = FALSE
	if (method=="Average") methodAverage = TRUE   # Required for later
	if (method!="function") methodFunction = NULL # Required for later

	if ( sum(rowGroup=="",na.rm=TRUE)>0 ){
	   warning(paste("rowGroup contains blanks. It is strongly recommended that you remove",
                    "these rows before calling the function.\n",
                    "   But for your convenience, the collapseRow function will remove these rows"));
	   rowGroup[rowGroup==""]=NA
	}

	# datET is a numeric matrix whose rows correspond to variables
	# e.g. probes of a microarray and whose columns to observations
	# e.g. microarrays 

	if ( sum(is.na(rowGroup))>0 ){
       warning(paste("The argument rowGroup contains missing data. It is strongly recommended\n",
              "   that you remove these rows before calling the function. Or redefine rowGroup\n",
			  "   so that it has no missing data. But for convenience, we remove these data."))
	}	
	
    ## Test to make sure the variables are the right length.
    #     if not, fix it if possible, or return 0 if not possible
	rowID  = as.character(rowID)
	rowGroup = as.character(rowGroup)
	rnDat = rownames(datET)
	if (length(rowID)!=length(rowGroup)){
		write("Error: rowGroup and rowID not the same length... exiting.","")
		return(0)
	}
	
	if (length(unique(rowID)) !=length(rowID) ){stop("rowID contains duplicate entries. Make sure that the argument rowID contains unique entries")}
	
	names(rowGroup) = rowID
	
    if ( sum(is.na(rowID))>0 ){warning("The argument rowID contains missing data. I recommend you choose non-missing, unique values for rowID, e.g. character strings.")}
	
	if ((is.null(rnDat))&(dim(datET)[1]==length(rowID))){
		write("Warning: *datET* does not have row names.  Assigning *rowID* as row names.","")
		rnDat <- rownames(datET) <- rowID
	}
	if (is.null(rnDat)){
		write("Error: *datET* does not have row names and length of *rowID*...","")
		write("... is not the same as # rows in *datET*... exiting.","")
		return(0)
	}
	if (sum(is.element(rnDat,rowID))!=length(rowID)){
		write("Warning: row names of input data and probes not identical...","")
		write("... Attempting to proceed anyway. Check results carefully.","")
		keepProbes = is.element(rowID, rownames(datET))
		rowID = rowID[keepProbes]
		datET= datET[rowID,]
		rowGroup = rowGroup[rowID]

	}

    restRows = (rowGroup!="" & !is.na(rowGroup))
    datET= datET[restRows,]
    rowGroup = rowGroup[restRows]
    rowID = rowID[restRows]
    rnDat = rnDat[restRows]


## For each group, select the row with the fewest missing values (if selectFewestMissing==TRUE)
##  Also, remove all rows with more than 90% missing data
	
	datET_in = datET  # This will be used as a reference later
		
## For each gene, select the gene with the fewest missing probes (if selectFewestMissing==TRUE)
##  Also, remove all probes with more than 90% missing data
	
	keep = .selectFewestMissing(datET, rowID, rowGroup, selectFewestMissing)
        datET = datET[keep, ];
        rowGroup = rowGroup[keep];
        rowID = rowID[keep];

	rnDat = rownames(datET)
	
##   If 0 < thresholdCombine < 1, only combine ids into their corresponding group if their 
#    correlation is greater than thresholdCombine.  This parameter supercedes all remaining 
#    parameters.
	
	if(!is.na(thresholdCombine)){
		if(!is.numeric(thresholdCombine)){
			write("thresholdCombine is not between -1 and 1 and is therefore being treated as NA","")
		} else if((thresholdCombine<(-1))|(thresholdCombine>1)){
			write("thresholdCombine is not between -1 and 1 and is therefore being treated as NA","")
		} else {
			output = .filterSimilarPS(datET, rowGroup, rowID, thresholdCombine)
			return(output)	
		}
	}
	
##   If method="function", use the function "methodFunction" as a way of combining genes
#    Alternatively, use one of the built-in functions 
#    Note: methodFunction must be a function that takes a vector of numbers as input and
#     outputs a single number. This function will return(0) or crash otherwise.

        recMethods = c("function","ME","MaxMean","maxRowVariance","MinMean","absMinMean","absMaxMean","Average");
        imethod = pmatch(method, recMethods);
        
	if (is.na(imethod)) {
		printFlush("Error: entered method is not a legal option. Recognized options are *maxRowVariance*,");
		printFlush("       *maxRowVariance*, *MaxMean*, *MinMean*, *absMaxMean*, *absMinMean*, *ME*,");
		printFlush("       *Average* or *function* for a user-defined function.")
		return(0)
	}
        if (imethod > 2) method = spaste(".", method);
	if (method=="function") 
        {
          method = methodFunction
	  if((!is.function(methodFunction))&(!is.null(methodFunction))){
            write("Error: *methodFunction* must be a function... please read the help file","")
            return(0)
	  }
        }
	if (!is.function(method)) if (method!="ME") method = get(method, mode = "function")
		
## Format the variables for use by this function
	rowID[is.na(rowID)] = rowGroup[is.na(rowID)]    # Use group if row is missing
	rownames(datET)[is.na(rnDat)]   = rowGroup[is.na(rnDat)]
	remove       = (is.na(rowID))|(is.na(rowGroup)) # Omit if both gene and probe are missing
	rowID  = rowID[!remove];
	rowGroup = rowGroup[!remove];
	names(rowGroup) = rowID
	rowID = sort(intersect(rnDat,rowID))
	if (length(rowID)<=1){
		write("Error: none of the *datET* rownames are in *rowID*...","")
		write("... please add rownames and try again... exiting.","")
		return(0)
	}
	rowGroup = rowGroup[rowID]
	datET  = as.matrix(datET)
	datET  = datET[rowID,]
	probes = rownames(datET)
	genes  = rowGroup[probes]
	tGenes = table(genes)
	datETOut=matrix(0,nrow=length(tGenes),ncol=ncol(datET))
	colnames(datETOut) = colnames(datET)
	rownames(datETOut) = sort(names(tGenes))
	rowsOut = rownames(datETOut)
	names(rowsOut) = rowsOut
	
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
	
	whichTestFn <- function(x){
		d    = datETOut[g,]
		test = (!is.na(x))&(!is.na(d))
		return(sum(x[test]==d[test]))
	}
	
# If method=ME, this function acts as the function moduleEigengene from the WGCNA library
	if (!is.function(method)) if (method=="ME"){
		datETOut = t(moduleEigengenes(t(datET),genes)$eigengenes)
		colnames(datETOut) = colnames(datET)
		rownames(datETOut) = substr(rownames(datETOut),3,nchar(rownames(datETOut)))
		out2 = cbind(rownames(datETOut),paste("ME",rownames(datETOut),sep="."))
		colnames(out2) = c("group","selectedRowID")
		out3 = is.element(rownames(datET_in),"@#$%^&*")
		names(out3) = rownames(datET_in)
		return(list(datETcollapsed = datETOut, group2row = out2, selectedRow = out3))		
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
		rowsOut[g] = probes[genes==g]
	}
        count = 0;
	for (g in twos){
		datETTmp = datET[probes[genes==g],]
		datETOut[g,] = as.numeric(method(datETTmp))
		whichTest    = apply(datETTmp,1,whichTestFn)
		rowsOut[g] = (names(whichTest)[whichTest==max(whichTest)])[1]
		count = count + 1;
		if (count %% 1000 == 0) collectGarbage();
	}
	for (g in more){
		datETTmp = datET[probes[genes==g],]
		adj = (0.5+0.5*cor(t(datETTmp),use="p"))^connectivityPower
		datETOut[g,] = as.numeric(datETTmp[which.max(rowSums(adj,na.rm=TRUE)),])
		whichTest    = apply(datETTmp,1,whichTestFn)
                rowsOut[g] = (names(whichTest)[whichTest==max(whichTest)])[1]
		count = count + 1;
		if (count %% 1000 == 0) collectGarbage();
	}
	if (!is.null(methodFunction))
		write("...Ignore previous comment.  Function completed properly!","")

		
# Retreive the information about which probes were saved, and include that information
#   as part of the output.  If method="function" or "Average" output placeholder values.
	if (!is.null(methodFunction)) {
		out2 = cbind(rownames(datETOut),paste("function",rownames(datETOut),sep="."))
		colnames(out2) = c("group","selectedRowID")
		out3 = is.element(rownames(datET_in),"@#$%^&*")
		names(out3) = rownames(datET_in)		
		return(list(datETcollapsed = datETOut, group2row = out2, selectedRow = out3))
	}
	if (methodAverage) {
		out2 = cbind(rownames(datETOut),paste("Average",rownames(datETOut),sep="."))
		colnames(out2) = c("group","selectedRowID")
		out3 = is.element(rownames(datET_in),"@#$%^&*")
		names(out3) = rownames(datET_in)		
		return(list(datETcollapsed = datETOut, group2row = out2, selectedRow = out3))
	}
	out2 = cbind(rownames(datETOut),rowsOut)
	colnames(out2) = c("group","selectedRowID")
	out3 = is.element(rownames(datET_in),rowsOut)
	names(out3) = rownames(datET_in)
	output = list(datETcollapsed = datETOut, group2row = out2, selectedRow = out3)
	return(output)
		
# End of function
} 
