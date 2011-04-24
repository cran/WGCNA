# Plots the kME values of genes in two groups of expression data for each module in an inputted color vector
kMEcomparisonScatterplot <- function (datExpr1, datExpr2, colorh, inA=NULL, inB=NULL, 
MEsA=NULL, MEsB=NULL, nameA="A", nameB="B", plotAll=FALSE, noGrey=TRUE, maxPlot=1000,
pch=19, 
fileName = if (plotAll) paste("kME_correlations_between_",nameA,"_and_",nameB,"_all.pdf",sep="") else
  paste("kME_correlations_between_",nameA,"_and_",nameB,"_inMod.pdf",sep=""),
...){

# First, get the data
	if (is.null(dim(datExpr1))) { 
		write ("Error: datExpr1 must be a matrix",""); return(0)
	}
	if (is.null(datExpr2)){
		datA = datExpr1[inA,]
		datB = datExpr1[inB,]
		if ((is.null(dim(datA)))|(is.null(dim(datB)))) { 
			 write ("Error: Check input for inA and inB.",""); return(0)
		}
	} else {
		if (is.null(dim(datExpr2))) { 
			write ("Error: datExpr2 must be a matrix",""); return(0)
		}
		datA = datExpr1
		datB = datExpr2
	}
	if ((dim(datA)[2]!=length(colorh))|(dim(datB)[2]!=length(colorh))){
		write ("Error: Both sets of input data and color vector must all have same length.",""); return(0)
	}
	
	if(is.null(MEsA))
		MEsA = (moduleEigengenes(datA, colors=as.character(colorh), excludeGrey=noGrey))$eigengenes
	if(is.null(MEsB))
		MEsB = (moduleEigengenes(datB, colors=as.character(colorh), excludeGrey=noGrey))$eigengenes
	mods  = substring(names(MEsA),3)
	kMEsA = as.data.frame(cor(datA,MEsA,use="p"))
	kMEsB = as.data.frame(cor(datB,MEsB,use="p"))
	
# Second, make the plots
	xlab  = paste("kME values in",nameA)
	ylab  = paste("kME values in",nameB)
        printFlush(paste("Plotting kME scatterplots into file", fileName));
	if (plotAll){
		pdf(file=fileName)
		numPlot = min(maxPlot,length(colorh));
		these   = sample(1:length(colorh),numPlot)
		for (i in 1:length(mods)){
			plotCol = mods[i]; if(mods[i]=="white") plotCol="black"
			verboseScatterplot(kMEsA[these,i],kMEsB[these,i],main=mods[i],
							   xlab=xlab,ylab=ylab,pch=pch,col=plotCol,...)
		}
		dev.off()
		return("DONE - Plotted All")
	}
	pdf(file=fileName)
	for (i in 1:length(mods)){
		these   = colorh==mods[i]
		plotCol = mods[i]; if(mods[i]=="white") plotCol="black"
		verboseScatterplot(kMEsA[these,i],kMEsB[these,i],main=mods[i],
						   xlab=xlab,ylab=ylab,pch=pch,col=plotCol,...)
	}
	dev.off()
	return("DONE - Plotted only in module")
}
