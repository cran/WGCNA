## Merge modules and reassign module genes based on kME values
moduleMergeUsingKME <- function (datExpr, colorh, ME=NULL, threshPercent=50, mergePercent = 25, reassignModules=TRUE, 
convertGrey=TRUE, omitColors="grey", reassignScale=1, threshNumber=NULL) {

## First assign all of the variables and put everything in the correct format.
	if (length(colorh)!=dim(datExpr)[2]){
		write("Error: color vector much match inputted datExpr columns.",""); return(0)
	}
	if (is.null(ME))  
		ME = (moduleEigengenes(datExpr, colors=as.character(colorh), excludeGrey=TRUE))$eigengenes
	if (dim(ME)[1]!=dim(datExpr)[1]){
		write("Error: ME rows much match inputted datExpr rows (samples).",""); return(0)
	}
	modules = colnames(ME)
	if (length(grep("ME",modules))==length(modules)) modules = substr(modules,3,nchar(modules))
	if (length(grep("PC",modules))==length(modules)) modules = substr(modules,3,nchar(modules))
	if (length(setdiff(modules,as.character(colorh)))>0){
		write("ME cannot include colors with no genes assigned.",""); return(0)
	}
	names(ME) = modules
	
	datCorrs = as.data.frame(cor(datExpr,ME,use="p"))
	colnames(datCorrs) = modules
	modules  = sort(modules[!is.element(modules,omitColors)])
	modulesI = modules # To test whether merging occurs
	datCorrs = datCorrs[,modules]
	iteration = 1
	if(is.null(colnames(datExpr))) colnames(datExpr) = as.character(1:length(colorh))
	rownames(datCorrs) <- colnames(datExpr)
	colorOut = colorh
	colorOut[is.element(colorOut,omitColors)] = "grey"
	mergeLog = NULL
	datExpr = t(datExpr) # For consistency with how the function was originally written.

## Iteratively run this function until no further changes need to be made
	while (!is.na(iteration)){
		write("",""); write("__________________________________________________","")
		write(paste("This is iteration #",iteration,". There are ",length(modules)," modules.",sep=""),"")
		iteration = iteration+1
		
	## Reassign modules if requested by reassignModules and convertGrey
		colorMax  = NULL
		whichMod  = apply(datCorrs,1,which.max)
		cutNumber = round(table(colorOut)*threshPercent/100)
		cutNumber = cutNumber[names(cutNumber)!="grey"] 
		if(!is.null(threshNumber)) cutNumber = rep(threshNumber,length(cutNumber))
		cutNumber = apply(cbind(cutNumber,10),1,max)
		for (i in 1:length(whichMod)) 
			colorMax = c(colorMax,modules[whichMod[i]])
		for (i in 1:length(modules)){
			corrs    = as.numeric(datCorrs[,i])
			cutValue = sort(corrs[colorOut==modules[i]],decreasing=TRUE)[cutNumber[i]]
			inModule = corrs>(cutValue*reassignScale)		
			if(convertGrey)
				colorOut[inModule&(colorOut=="grey")&(colorMax==modules[i])]=modules[i]
			if(reassignModules)
				colorOut[inModule&(colorOut!="grey")&(colorMax==modules[i])]=modules[i]
		}
		
	## Merge all modules meeting the mergePercent and threshPercent criteria
		for (i in 1:length(modules)){
			cutNumber  = round(table(colorOut)*threshPercent/100)
			cutNumber  = cutNumber[names(cutNumber)!="grey"]
			if(!is.null(threshNumber)) cutNumber = rep(threshNumber,length(cutNumber))
			cutNumber = apply(cbind(cutNumber,10),1,max)
			corrs      = as.numeric(datCorrs[,i])
			# Make sure you do not include more genes than are in the module
			numInMod   = sum(colorOut==modules[i])
			cutValue   = sort(corrs[colorOut==modules[i]],decreasing=TRUE)[min(numInMod,cutNumber[modules[i]])]
			colorMod   = colorOut[corrs>=cutValue]
			colorMod   = colorMod[colorMod!="grey"]
			modPercent = 100*table(colorMod)/length(colorMod)
			modPercent = modPercent[names(modPercent)!=modules[i]]
			if(length(modPercent)>1) if(max(modPercent)>mergePercent){
				whichModuleMerge = names(modPercent)[which.max(modPercent)]
				colorOut[colorOut==modules[i]] = whichModuleMerge
				write(paste(modules[i],"has been merged into",whichModuleMerge,"."),"")
				mergeLog = rbind(mergeLog,c(modules[i],whichModuleMerge))
			}
		}
		
	## If no modules were merged, then set iteration to NA 
		modules  = sort(unique(colorOut))		
		modules  = modules[modules!="grey"]
		if (length(modules)==length(modulesI)) iteration=NA
		modulesI = modules
		
	## Recalculate the new module membership values
		MEs      = (moduleEigengenes(t(datExpr), colors=as.character(colorOut)))$eigengenes
		MEs      = MEs[,colnames(MEs)!="MEgrey"]
		datCorrs = as.data.frame(cor(t(datExpr),MEs,use="p"));
		colnames(datCorrs) = modules
	}
	if(!is.null(dim(mergeLog))){
		colnames(mergeLog) = c("Old Module","Merged into New Module")
		rownames(mergeLog) = paste("Merge #",1:dim(mergeLog)[1],sep="")
	}
	return(list(moduleColors=colorOut,mergeLog=mergeLog))
}

