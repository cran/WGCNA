# Measure enrichment between inputted lists of genes (ie, modules) and user-defined lists.
userListEnrichment <- function (geneR, labelR, fnIn=NULL, catNmIn=fnIn, 
nameOut = "enrichment.csv",useBrainLists=FALSE, omitCategories ="grey",
outputCorrectedPvalues=TRUE){

# Quickly test for input compatability
	if (length(geneR)!=length(labelR)){
		write("ERROR: geneR and labelR must have same number of elements.",""); return(0)
	}
	if (length(catNmIn)<length(fnIn)){
		catNmIn=c(catNmIn,fnIn[(length(catNmIn)+1):length(fnIn)])
		write("WARNING: not enough category names.  
			   Naming remaining categories with file names.","")
	}
	if (is.null(fnIn)&(!useBrainLists)){
		write("ERROR: Either enter user-defined lists or set useBrainLists=TRUE.","")
		return(0)
	}
	
# Read in the data
	glIn = NULL
	for (i in 1:length(fnIn)){
		ext = substr(fnIn[i],nchar(fnIn[i])-2,nchar(fnIn[i]))
		if (ext=="csv"){
			datIn = read.csv(fnIn[i])
			if (colnames(datIn)[2]=="Gene") {
				datIn = datIn[,2:3]
			} else { datIn = datIn[,1:2] }
		} else{
			datIn = scan(fnIn[i], what="character",sep="\n")
			datIn = cbind(datIn[2:length(datIn)],datIn[1])
		}
		colnames(datIn)=c("Gene","Category")
		datIn[,2] = paste(datIn[,2],catNmIn[i],sep="__")
		glIn = rbind(glIn,datIn)
	}
	if(useBrainLists) { 
		# load("BrainLists.RData") # *** Set path to correct folder ***
                BrainLists = 0; # This is to satisfy R CMD check complaining about no visible binding
                data(BrainLists)
		write("See help file for details regarding brain list references.","")
		glIn = rbind(glIn, BrainLists)
	}
	# Remove any duplicate entries
	removeDups = unique(paste(as.character(glIn[,1]),as.character(glIn[,2]),sep="@#$%"))
	if (length(removeDups)<length(glIn[,1]))
		glIn = t(as.matrix(as.data.frame(strsplit(removeDups,"@#$%",fixed=TRUE))))
	geneIn  = as.character(glIn[,1])
	labelIn = as.character(glIn[,2])
	
# Format the data, excluding anything not in the inputted gene list 
	geneAll = sort(unique(geneR))
	keep    = is.element(geneIn,geneAll)
	geneIn  = geneIn[keep]
	labelIn = labelIn[keep]
	catsR   = sort(unique(labelR))
	catsR   = catsR[!is.element(catsR,omitCategories)]
	catsIn  = sort(unique(labelIn))
	lenAll  = length(geneAll)
	
# Determine the hypergeometric enrichment values 
	results = list(pValues = NULL, ovGenes = list(), sigOverlaps=NULL)
	namesOv = NULL
	for(r in 1:length(catsR)) for (i in 1:length(catsIn)) {
		gr = geneR[labelR==catsR[r]]
		gi = geneIn[labelIn==catsIn[i]]
		go = intersect(gr,gi)
		lr = length(gr)
		li = length(gi)
		lo = length(go)
		pv = .phyper2(lenAll,lr,li,lo,FALSE)
		if (lo==0) pv=1
		if (pv<0.0001) 	pv = .phyper2(lenAll,lr,li,lo,TRUE)
		pOut = c(catsR[r],catsIn[i],lo,pv)
		results$pValues = rbind(results$pValues,pOut)
		namesOv = c(namesOv,paste(catsR[r],"--",catsIn[i]))
		results$ovGenes[[length(namesOv)]] = go
	}
	results$pValues = cbind(results$pValues, 
		apply(cbind(1,as.numeric(results$pValues[,4])*length(namesOv)),1,min))
	colnames(results$pValues) = c("InputCategories","UserDefinedCategories",
								  "NumOverlap","Pvalues","CorrectedPvalues")
	names(results$ovGenes) = namesOv
	results$sigOverlaps = results$pValues[as.numeric(results$pValues[,5])<0.05,c(1,2,5)]
	if(!outputCorrectedPvalues){
		results$sigOverlaps = results$pValues[as.numeric(results$pValues[,4])<0.05,c(1,2,4)]
		write("Note that outputted p-values are not corrected for multiple comparisons.","")
	}
	results$sigOverlaps = results$sigOverlaps[order(as.numeric(results$sigOverlaps[,3])),]
	
	write.csv(results$sigOverlaps,nameOut, row.names=FALSE)
	write(paste(length(namesOv),"comparisons were successfully performed."),"")
	return(results)
}

.phyper2 <- function (total, group1, group2, overlap, verySig=TRUE ,lt=TRUE){
# This function is the same is phyper, just allows for more sensible input values
	q = overlap
	m = group1
	n = total-group1
	k = group2
	prob = phyper(q, m, n, k, log.p = verySig, lower.tail=lt)
	if (verySig) return(-prob)
	return(1-prob)
}
