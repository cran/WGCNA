# Makes a consensus network using all of the default values in the WGCNA library.

consensusDissTOMandTree <- function (multiExpr, softPower, TOM=NULL){
	nGenes = dim(multiExpr[[1]]$data)[2]
	nSets  = length(multiExpr)
	if(is.null(TOM)){
		adjacencies <- TOM <- list()
		for (set in 1:nSets){
			adjacencies[[set]] = adjacency(multiExpr[[set]]$data,power=softPower,type="signed");
			diag(adjacencies[[set]])=0
			write(paste("Adjacency, set",set),"")
			TOM[[set]] = TOMsimilarity(adjacencies[[set]], TOMType="signed");
			write(paste("Similarity, set",set),"")
			collectGarbage()
		}
	}
	nSets  = length(TOM)
	set.seed(12345);     
	scaleP = 0.95;
	nSamples = as.integer(1/(1-scaleP) * 1000);
	scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
	TOMScalingSamples = list();
	scaleQuant <- scalePowers <- rep(1, nSets)
	for (set in 1:nSets){
		TOMScalingSamples[[set]] = as.dist(TOM[[set]])[scaleSample]
		scaleQuant[set] = quantile(TOMScalingSamples[[set]],probs = scaleP, type = 8);
		if (set>1){
			scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
			TOM[[set]] = TOM[[set]]^scalePowers[set];
		}
		write(paste("Scaling, set",set),"")
	}
	
	half = round(nGenes/2); haP1 = half+1
	kp   = list(list(c(1:half),c(1:half)),list(c(1:half),c(haP1:nGenes)),
		  	    list(c(haP1:nGenes),c(1:half)),list(c(haP1:nGenes),c(haP1:nGenes)))
	consensusTOMi = list()
	for (i in 1:4){
		a = kp[[i]][[1]];  b = kp[[i]][[2]]
		consensusTOMi[[i]] = TOM[[1]][a,b]
		for (j in 2:nSets)   consensusTOMi[[i]] = pmin(consensusTOMi[[i]], TOM[[j]][a,b]);
		write(paste(i,"of 4 iterations in pMin"),"")
	}
	consensusTOM = rbind(cbind(consensusTOMi[[1]],consensusTOMi[[2]]), 
						 cbind(consensusTOMi[[3]],consensusTOMi[[4]]))
	rownames(consensusTOM) <- colnames(consensusTOM) <- colnames(multiExpr[[1]]$data)
	
	consensusTOM = 1-consensusTOM
	write("Starting dendrogram tree.","")
	consTree     = fastcluster::hclust(as.dist(consensusTOM), method = "average");
	write("DONE!!!!","")
	out = list(consensusTOM,consTree)
	names(out) = c("consensusTOM","consTree")
	return(out)
}

.collect_garbage <- function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}
