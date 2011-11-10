populationMeansInAdmixture<-function (
datProportions, 
datE.Admixture, 
scaleProportionsTo1 = TRUE, 
scaleProportionsInCelltype=TRUE,
setMissingProportionsToZero = FALSE
) 
{
    datProportions = as.matrix(datProportions)
    if (dim(datE.Admixture)[[1]] != dim(as.matrix(datProportions))[[1]]) 
        stop("Input error. The numbers of samples are not congruent: dim(datE.Admixture)[[1]] is unequal to dim(datProportions)[[1]]. Hint: Consider transposing one of the matrices.")
    noMissing = apply(is.na(datProportions), 1, sum)
    if (max(noMissing) > 0) {
        warning(paste("Urgent Warning: datProportions contains missing proportions in the following row(s):", 
            paste(which(noMissing > 0), collapse = ","), "\nCheck these rows in datProportions. But for your convenience, I will proceed"))
    }
    if (setMissingProportionsToZero) {
        datProportions[is.na(datProportions)] = 0
    }
    noNegative = apply((datProportions) < 0, 1, sum, na.rm = T)
    if (max(noNegative) > 0) {
        stop(paste("datProportions contains negative numbers. Negative proportions can be found in the following row(s):", 
            paste(which(noNegative > 0), collapse = ","), "\nCheck these rows in datProportions."))
    }
    sumsTo1 = TRUE
    sumProp = apply(datProportions, 1, sum, na.rm = TRUE)
    if (max(sumProp, na.rm = T) > 1.0000001) {
        sumsTo1 = FALSE
        if (scaleProportionsTo1) {
            warning(paste("The sum of proportions is larger than 1 for some rows including:", 
                paste(which(sumProp > 1.0000001)[1:5], collapse = ","), 
                "\nCheck these rows in datProportions. By default, I will scale them so that they sum to 1.\nBut if you do not want this scaling, please set scaleProportionsTo1=FALSE .\n"))
        }
    }
    if (min(sumProp, na.rm = T) < 0.99999) {
        sumsTo1 = FALSE
        if (scaleProportionsTo1) {
            warning(paste("The sum of proportions is smaller than 1 for some rows including:", 
                paste(which(sumProp < 0.99999)[1:5], collapse = ","), 
                "\n Check these rows in datProportions. By default, I will scale them so that they sum to 1.\nBut if you do not want this scaling, please set scaleProportionsTo1=FALSE .\n"))
        }
    }
    if (scaleProportionsTo1) {
        sumsTo1 = TRUE
        for (i in 1:dim(datProportions)[1]) {
            datProportions[i, ] = datProportions[i, ]/sum(datProportions[i, 
                ], na.rm = T)
        }
    }
    if (sumsTo1) {
	if(scaleProportionsInCelltype) {
		for(ci in dim(datProportions)[2]) 
		datProportions[,ci]=datProportions[,ci]-mean(datProportions[,ci])	
	}

        fit1 = lm(as.matrix(datE.Admixture) ~ . - 1, data = data.frame(datProportions))
        datPredictedMeans = t(as.matrix(fit1$coefficients))
    }
    if (!sumsTo1) {
	if(scaleProportionsInCelltype) {
		for(ci in dim(datProportions)[2]) 
		datProportions[,ci]=datProportions[,ci]-mean(datProportions[,ci])	
	}
        fit1 = lm(as.matrix(datE.Admixture) ~ ., data = data.frame(datProportions))
        if (dim(as.matrix(datProportions))[[2]] == 1) {
            datPredictedMeans = (matrix(fit1$coefficients[-1, 
                ], ncol = 1))
        }
        else {
            datPredictedMeans = t(as.matrix(fit1$coefficients[-1, 
                ]))
        }
    }
    dimnames(datPredictedMeans)[[1]] = dimnames(datE.Admixture)[[2]]
    if (is.null(dimnames(datPredictedMeans)[[2]])) {
        dimnames(datPredictedMeans)[[2]] = paste("Mean", 1:dim(datPredictedMeans)[[2]], 
            sep = ".")
    }
    else {
        dimnames(datPredictedMeans)[[2]] = paste("Mean", dimnames(datPredictedMeans)[[2]], 
            sep = ".")
    }
    datPredictedMeans
}
