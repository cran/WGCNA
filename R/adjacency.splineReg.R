adjacency.splineReg = function(datExpr, df = 6-(nrow(datExpr)<100)-(nrow(datExpr)<30), symmetrizationMethod = "mean", ...) {
  
  if (!is.element(symmetrizationMethod, c("none", "min" ,"max", "mean"))) {
    stop("Unrecognized symmetrization method.")
  }

  datExpr = matrix(as.numeric(as.matrix(datExpr)), nrow(datExpr), ncol(datExpr))
  n = ncol(datExpr)
  splineRsquare = matrix(NA, n,n)

  for (i in 2:n) {
    for (j in 1:(i-1)) {
      del = is.na(datExpr[, i]+datExpr[,j]) 
      if (sum(del)>=(n-1) | var(datExpr[, i], na.rm=T)==0 | var(datExpr[, j], na.rm=T)==0) {
	  splineRsquare[i, j] = splineRsquare[j, i]=NA
	}else{
	  dati = datExpr[!del, i]; datj = datExpr[!del, j];
        lmSij=glm( dati ~ ns( datj, df = df, ...))
        splineRsquare[i, j] = cor( dati, predict(lmSij))^2
        lmSji=glm( datj ~ ns( dati, df = df, ...))
        splineRsquare[j, i] = cor( datj, predict(lmSji))^2
        rm(dati, datj, lmSij, lmSji)
      }
    }
  }

  diag(splineRsquare) = rep(1,n)
  if (symmetrizationMethod =="none") {adj= splineRsquare} else
 {    adj = switch(symmetrizationMethod, 
	min = pmin(splineRsquare, t(splineRsquare)),
	max = pmax(splineRsquare, t(splineRsquare)),
	mean = (splineRsquare + t(splineRsquare))/2)}

  adj

}


