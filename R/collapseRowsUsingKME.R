# This function chooses a single plobe per gene based on a kME table
collapseRowsUsingKME <- function (MM, Gin, Pin=NULL, kMEcols = 1:dim(MM)[2]){
	
	if (is.null(Pin))  Pin = rownames(MM)
	rownames(MM) = Pin
	Gout = as.character(sort(unique(Gin)))
	cors = MM[,kMEcols]
	maxC = apply(cors,1,max)
	
	MMout = matrix(0, nrow = length(Gout), ncol = dim(MM)[2])
	colnames(MMout) = colnames(MM)
	rownames(MMout) = Gout
	MM = as.matrix(MM)
	
	keepThese = NULL
	for (g in 1:length(Gout)){
		maxCg = maxC
		maxCg[Gin!=Gout[g]] = -1000
		keep      = which(maxCg==max(maxCg))[1]
		MMout[g,] = MM[keep,]
		keepThese = c(keepThese, keep)
	}
	group2Row   = cbind(Gout,Pin[keepThese])
	colnames(group2Row) = c("group","selectedRowID")
	selectedRow = is.element(1:length(Pin),keepThese)
	out = list(MMout, group2Row, selectedRow)
	names(out) = c("MMcollapsed","group2Row", "selectedRow")
	return(out)
}

