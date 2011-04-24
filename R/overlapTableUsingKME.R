# Determines significant overlap between modules in two networks based on kME tables.
overlapTableUsingKME <- function(dat1, dat2, colorh1, colorh2, MEs1=NULL, MEs2=NULL, 
name1="MM1", name2="MM2", cutoffMethod="assigned", cutoff=0.5, omitGrey=TRUE, 
datIsExpression=TRUE){
	
# Run a few tests on the imput data formatting
	if (is.null(dim(dat1))|is.null(dim(dat2))) { 
		write("Error: dat1 and dat2 must be matrices.",""); return(0)
	}
	if ((dim(dat1)[datIsExpression+1]!=length(colorh1))|
		(dim(dat2)[datIsExpression+1]!=length(colorh2))){
		write("Error: Both sets of input data and color vectors must have same length.",""); return(0)
	}
	if ((cutoffMethod=="pvalue")&(datIsExpression==FALSE)){
		write("Error: Pvalues are not calculated if datIsExpression=TRUE.  Choose other cutoffMethod.",
			  ""); return(0)
	}
	
# Find and format the kME values and other variables for both inputs
	G1 = dimnames(dat1)[[datIsExpression+1]];  G2 = dimnames(dat2)[[datIsExpression+1]];
	if(datIsExpression){
		if(is.null(MEs1))
			MEs1 = (moduleEigengenes(dat1, colors=as.character(colorh1), excludeGrey=omitGrey))$eigengenes
		if(is.null(MEs2))
			MEs2 = (moduleEigengenes(dat2, colors=as.character(colorh2), excludeGrey=omitGrey))$eigengenes
		mods1 = colnames(MEs1);  mods2 = colnames(MEs2)
		if (length(grep("ME",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("PC",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("ME",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))	
		if (length(grep("PC",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))
		out = corAndPvalue(dat1,MEs1);  MM1 = out$cor;  PV1 = out$p;  rm(out);
		out = corAndPvalue(dat2,MEs2);  MM2 = out$cor;  PV2 = out$p;  rm(out);
		colnames(MM1) <- colnames(PV1) <- mods1;
		colnames(MM2) <- colnames(PV2) <- mods2;
		rownames(MM1) <- rownames(PV1) <- G1;
		rownames(MM2) <- rownames(PV2) <- G2;
	} else {
		MM1 = dat1[,sort(colnames(dat1))];  mods1 = colnames(MM1)
		MM2 = dat2[,sort(colnames(dat2))];  mods2 = colnames(MM2)
		if (length(grep("ME",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("PC",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("ME",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))
		if (length(grep("PC",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))
		colnames(MM1) = mods1;  colnames(MM2) = mods2;
		rownames(MM1) = G1;     rownames(MM2) = G2;
		if(omitGrey){
			MM1 = MM1[,!is.element(mods1,"grey")];  mods1 = colnames(MM1)
			MM2 = MM2[,!is.element(mods2,"grey")];  mods2 = colnames(MM2)
		}
	}
	if ((length(setdiff(mods1,as.character(colorh1)))>omitGrey)|
		(length(setdiff(mods2,as.character(colorh2)))>omitGrey)){
		write("MEs cannot include colors with no genes assigned.",""); return(0)
	}
	l1 = length(mods1);	 l2 = length(mods2)			
	cutoffMethod = substr(cutoffMethod,1,1)
	names=c(name1,name2)
	comGenes = sort(unique(intersect(G1,G2)));  total   = length(comGenes)
	MM1     = MM1[comGenes,];  MM2 = MM2[comGenes,]
    if (datIsExpression){
		PV1 = PV1[comGenes,];  PV2 = PV2[comGenes,]
	}
	names(colorh1) = G1;  colorh1 = colorh1[comGenes]
	names(colorh2) = G2;  colorh2 = colorh2[comGenes]
	
# Assign each gene in each module to a vector corresponding to the modules
	genes1 <- genes2 <- list()
	if (cutoffMethod=="a"){
		for (i in 1:l1)  genes1[[i]] = comGenes[colorh1==mods1[i]]
		for (i in 1:l2)  genes2[[i]] = comGenes[colorh2==mods2[i]]
	} else if (cutoffMethod=="p") {
		for (i in 1:l1)  genes1[[i]] = comGenes[PV1[,mods1[i]]<=cutoff]
		for (i in 1:l2)  genes2[[i]] = comGenes[PV2[,mods2[i]]<=cutoff]
	} else if (cutoffMethod=="k") {
		for (i in 1:l1)  genes1[[i]] = comGenes[MM1[,mods1[i]]>=cutoff]
		for (i in 1:l2)  genes2[[i]] = comGenes[MM2[,mods2[i]]>=cutoff]
	} else if (cutoffMethod=="n") {
		for (i in 1:l1)  genes1[[i]] = comGenes[rank(-MM1[,mods1[i]])<=cutoff]
		for (i in 1:l2)  genes2[[i]] = comGenes[rank(-MM2[,mods2[i]])<=cutoff]
	} else {
		write("ERROR: cutoffMethod entered is not supported.",""); return(0)
	}
	names(genes1) = paste(names[1],mods1,sep="_")
	names(genes2) = paste(names[2],mods2,sep="_")
	
# Determine signficance of each comparison and write out all of the gene lists
	ovGenes = list()
	ovNames = rep("",l1*l2)
	pVals   = matrix(1, nrow=l1, ncol=l2)
	rownames(pVals) = paste(names[1],mods1,sep="_")
	colnames(pVals) = paste(names[2],mods2,sep="_")
	i = 0
	for (m1 in 1:l1) for (m2 in 1:l2) {
		i = i+1
		ovGenes[[i]] = sort(unique(intersect(genes1[[m1]],genes2[[m2]])))
		pVals[m1,m2] = .phyper2(total,length(genes1[[m1]]), length(genes2[[m2]]),length(ovGenes[[i]]))
		if (pVals[m1,m2]>10^(-10))   pVals[m1,m2] = 
		.phyper2(total,length(genes1[[m1]]), length(genes2[[m2]]),length(ovGenes[[i]]),FALSE)
		ovNames[i] = paste(names[1],mods1[m1],names[2],mods2[m2],sep="_")
	}
	names(ovGenes) = ovNames
	out = list(pVals,comGenes,genes1,genes2,ovGenes)
	names(out) = c("PvaluesHypergeo","AllCommonGenes",paste("Genes",names,sep=""),"OverlappingGenes")
	return(out)
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