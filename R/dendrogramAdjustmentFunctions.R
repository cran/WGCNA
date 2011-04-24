## This file contains several functions which can be used to adjust the dendrogram
#   in ways which keep the dendrogram mathematically identical (ie, branch swapping,
#   branch reflection, etc).  The goal is to biologically optimize the dendrogram.

# ----------------------------------------------------------------------------- #

orderBranchesUsingHubGenes <- function(hierTOM, datExpr=NULL, colorh=NULL, type="signed", 
adj=NULL, iter=NULL, useReflections=FALSE, allowNonoptimalSwaps=FALSE){

# First read in and format all of the variables
	hGenes = hierTOM$labels
	if(is.null(adj)){
		genes = chooseOneHubInEachModule(datExpr, colorh, type=type)
		adj   = adjacency(datExpr[,genes],type=type, power=(2-(type=="unsigned")))
		colnames(adj) <- rownames(adj) <- genes
	}
	genes = rownames(adj)
	if(length(genes)!=length(intersect(genes,hGenes))){
		write("All genes in the adjacency must also be in the gene tree.  Check to make sure","")
		write("that names(hierTOM$labels) is set to the proper gene or probe names and that","")
		write("these correspond to the expression / adjacency gene names.","")
		return(0)
	}
	genes = hGenes[is.element(hGenes,genes)]
	adj   = adj[genes,genes]
	if (is.null(iter)) iter = length(genes)^2
	iters=(1:iter)/iter
	swapAnyway = rep(0,length(iters))  # Quickly decreasing chance of random swap
	if (allowNonoptimalSwaps)  swapAnyway = ((1-iters)^3)/3+0.001
		
# Iterate random swaps in the branch, only accepting the new result
#  if it produces a higher correlation than the old result OR if the 
#  random variable says to swap (which gets less likely each iteration)
	changes=NULL
	for (i in 1:iter){
		swap = 1; 
		if (useReflections) swap = sample(0:1,1)
		gInd = sample(1:length(genes),2)
		g    = genes[gInd]
		if (swap==1) {
			hierTOMnew = swapTwoBranches(hierTOM, g[1], g[2])
		} else hierTOMnew = reflectBranch(hierTOM, g[1], g[2], TRUE)
		oldSum    = .offDiagonalMatrixSum(adj)
		oGenesNew = hGenes[hierTOMnew$order]
		oGenesNew = oGenesNew[oGenesNew%in%genes]
		adjNew    = adj[oGenesNew,oGenesNew]
		newSum    = .offDiagonalMatrixSum(adjNew)
		if ((newSum>oldSum)|((sample(1:1000,1)/1000)<swapAnyway[i])) {
			hierTOM = hierTOMnew
			changes = rbind(changes,c(i,ifelse(swap==1,"Swap","Reflect"),g,oldSum,newSum))
			adj     = adjNew
		}
		write(paste("Interation",i,"of",iter),"")
		collectGarbage()
	}
	
# Perform all of the suggested swappings on the input network.
	
# Output the results
	colnames(changes)=c("Iter.#","Swap?","Gene1","Gene2","OldScore","NewScore")
    out = list(geneTree = hierTOM, changeLog = changes)
	return(out)
}

# ----------------------------------------------------------------------------- #

selectBranch <- function (hierTOM, g1, g2){
## This function selects of all genes in a branch given a gene in the
##  branch (g1) and a gene in a neighboring branch (g2), returning the 
##  indices for genes in the branch in the hierTOM$labels vector 
	
# Convert genes to UNORDERED indices (if given, indices should be ordered)
	if(is.numeric(g1)) g1 = hierTOM$order[g1]
	if(is.numeric(g2)) g2 = hierTOM$order[g2]
	if(!is.numeric(g1)) g1 = which(hierTOM$labels==g1)
	if(!is.numeric(g2)) g2 = which(hierTOM$labels==g2)
	if((length(g1)==0)|(length(g2)==0)|(max(c(g1,g2))>length(hierTOM$labels))){
		write("Input genes are not both legal indices","")
		return(hierTOM);
	}
	
# Now determine which branch is the correct one, and find the genes
	len = length(hierTOM$height)
	tree1 = which(hierTOM$merge==(-g1))%%len
	continue=length(which(hierTOM$merge==tree1))>0
	while(continue){
		nextInd = which(hierTOM$merge==tree1[length(tree1)])%%len
		tree1 = c(tree1,nextInd)
		continue=length(which(hierTOM$merge==nextInd))>0
	}
	
	branchIndex = which(hierTOM$height==.minTreeHeight(hierTOM,g1,g2))
	branch=hierTOM$merge[branchIndex,]
	b1 <- NULL
	if(is.element(branch[1],tree1)){
		b1 = .getBranchMembers(hierTOM,branch[1],b1)
	} else b1 = .getBranchMembers(hierTOM,branch[2],b1)
	collectGarbage()
	return(b1)
}

# ----------------------------------------------------------------------------- #

reflectBranch <- function (hierTOM, g1, g2, both=FALSE){
## This function reverses the ordering of all genes in a branch of the
##  clustering tree defined by the minimal branch possible that contains
##  both g1 and g2 (as either ORDERED index or gene names), or just by
##  the genes in g1
	
	b1 = selectBranch(hierTOM, g1, g2)
	if (both) b1 = c(b1,selectBranch(hierTOM, g2, g1))
	
# Now reorder the hierTOM correctly
	ord = hierTOM$order
	i1 = which(ord%in%b1)
	b=1:(min(i1)-1); 
	if(b[length(b)]<b[1]) b = NULL
	e=(max(i1)+1):length(ord); 
	if((max(i1)+1)>length(ord)) e = NULL
	ord = ord[c(b,i1[order(i1,decreasing=T)],e)]
	hierTOM$order = ord
	return(hierTOM)
}

# ----------------------------------------------------------------------------- #

swapTwoBranches <- function (hierTOM, g1, g2){
## This function re-arranges two branches in a heirarchical clustering tree
##  at the nearest branch point of two given genes (or indices) 
	
# Convert genes to indices (ORDERED AS ON THE PLOT)
	if(is.numeric(g1)) g1 = hierTOM$order[g1]
	if(is.numeric(g2)) g2 = hierTOM$order[g2]
	if(!is.numeric(g1)) g1 = which(hierTOM$labels==g1)
	if(!is.numeric(g2)) g2 = which(hierTOM$labels==g2)
	if((length(g1)==0)|(length(g2)==0)|(max(c(g1,g2))>length(hierTOM$labels))){
		write("Input genes are not both legal indices","")
		return(hierTOM);
	}
	
# Now determine the genes in each branch
	branchIndex = which(hierTOM$height==.minTreeHeight(hierTOM,g1,g2))
	b1 <- b2 <- NULL
	b1 = .getBranchMembers(hierTOM,hierTOM$merge[branchIndex,1],b1)
	b2 = .getBranchMembers(hierTOM,hierTOM$merge[branchIndex,2],b2)
	
# Now reorder the hierTOM correctly
	ord = hierTOM$order
	i1 = which(ord%in%b1)
	i2 = which(ord%in%b2)
	if(min(i1)>min(i2)) {tmp = i1; i1=i2; i2=tmp; rm(tmp)}
	b=1:(min(i1)-1); 
	if(b[length(b)]<b[1]) b = NULL
	e=(max(i2)+1):length(ord); 
	if((max(i2)+1)>length(ord)) e = NULL
	ord = ord[c(b,i2,i1,e)]
	hierTOM$order = ord
	return(hierTOM)
}

# ----------------------------------------------------------------------------- #

chooseOneHubInEachModule <- function(datExpr, colorh, numGenes=100, 
omitColors="grey", power=2, type="signed",...){
## This function returns the gene in each module with the highest connectivity, given
#   a number of randomly selected genes to test.
	
	numGenes = max(round(numGenes),2)
	keep     = NULL
	isIndex  = FALSE
	modules  = names(table(colorh));
	numCols  = table(colorh)
	if(!(is.na(omitColors)[1]))  modules = modules[!is.element(modules,omitColors)]
	if(is.null(colnames(datExpr))){
		colnames(datExpr) = 1:dim(datExpr)[2]
		isIndex = TRUE
	}
	
	for (m in modules){
		num   = min(numGenes,numCols[m])
		inMod = which(is.element(colorh,m)) 
		keep  = c(keep, sample(inMod,num))
	}
	colorh  = colorh[keep]
	datExpr = datExpr[,keep]
	return(chooseTopHubInEachModule(datExpr, colorh, omitColors, power, type,...))
}

# ----------------------------------------------------------------------------- #

chooseTopHubInEachModule <- function(datExpr, colorh, omitColors="grey", 
power=2, type="signed",...){
## This function returns the gene in each module with the highest connectivity.
	
	isIndex = FALSE
	modules = names(table(colorh));
	if(!(is.na(omitColors)[1]))  modules = modules[!is.element(modules,omitColors)]
	if(is.null(colnames(datExpr))){
		colnames(datExpr) = 1:dim(datExpr)[2]
		isIndex = TRUE
	}
	
	hubs = rep(NA,length(modules))
	names(hubs) = modules
	for (m in modules){
		adj = adjacency(datExpr[,colorh==m],power=power,type=type,...)
		hub = which.max(rowSums(adj))
		hubs[m] = colnames(adj)[hub]
	}
	if (isIndex){
		hubs = as.numeric(hubs)
		names(hubs) = modules
	}
	return(hubs)
}

#################################################################################
# Internal functions.............................................................

options(expressions=50000) # Required for .getBranchMembers

.getBranchMembers <- function(hierTOM, ind, members){
# This is a recursive function that gets all the indices of members of
#  a branch in an hClust tree.
	if(ind<0) return(c(members,-ind))
	m1 = hierTOM$merge[ind,1]
	m2 = hierTOM$merge[ind,2]
	if (m1>0) {
		members = .getBranchMembers(hierTOM,m1,members)
	} else members = c(members,-m1)
	if (m2>0) {
		members = .getBranchMembers(hierTOM,m2,members)
	} else members = c(members,-m2)
	return(members)
}

# ----------------------------------------------------------------------------- #

.minTreeHeight <- function(hierTOM,l1,l2) {
## This function finds the minimum height at which two leafs
##  in a hierarchical clustering tree are connected.  l1 and
##  l2 are the UNORDERED indices for the two leafs.
	
## Return 2 (larger than 1, if l1 or l2 is negative). This represents 
##  positions that are off the edge of the tree.
	if((l1<0)|(l2<0)) return(2)
	
## Get the tree for l1
	len = length(hierTOM$height)
	tree1 = which(hierTOM$merge==(-l1))%%len
	continue=length(which(hierTOM$merge==tree1))>0
	while(continue){
		nextInd = which(hierTOM$merge==tree1[length(tree1)])%%len
		tree1 = c(tree1,nextInd)
		continue=length(which(hierTOM$merge==nextInd))>0
	}
	
## Get the tree for l2
	tree2 = which(hierTOM$merge==(-l2))%%len
	continue=length(which(hierTOM$merge==tree2))>0
	while(continue){
		nextInd = which(hierTOM$merge==tree2[length(tree2)])%%len
		tree2 = c(tree2,nextInd)
		continue=length(which(hierTOM$merge==nextInd))>0
	}
	
## Now find the index where the two trees first agree
	minTreeLen = min(c(length(tree1),length(tree2)))
	tree1 = tree1[(length(tree1)-minTreeLen+1):length(tree1)]
	tree2 = tree2[(length(tree2)-minTreeLen+1):length(tree2)]
	treeInd = tree1[min(which(tree1==tree2))]
	
## Now find and return the minimum tree height
	return(hierTOM$height[ifelse(treeInd==0,len,treeInd)])
}

# ----------------------------------------------------------------------------- #

.offDiagonalMatrixSum <- function(adj){
    len = dim(adj)[1]	
	output=sum(diag(adj[1:(len-1),2:len]))
	return(output)	
}
