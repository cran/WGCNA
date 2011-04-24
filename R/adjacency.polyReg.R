adjacency.polyReg = function(datExpr, degree=3, symmetrizationMethod = "mean") {
  
  if (!is.element(symmetrizationMethod, c("none", "min" ,"max", "mean"))) {
    stop("Unrecognized symmetrization method.")
  }

  datExpr = matrix(as.numeric(as.matrix(datExpr)), nrow(datExpr), ncol(datExpr))
  n = ncol(datExpr)
  polyRsquare = matrix(NA, n,n)

  for (i in 2:n) {
    for (j in 1:(i-1)) {
      del = is.na(datExpr[, i]+datExpr[,j]) 
      if (sum(del)>=(n-1) | var(datExpr[, i], na.rm=T)==0 | var(datExpr[, j], na.rm=T)==0) {
	  polyRsquare[i, j] = polyRsquare[j, i]=NA
	}else{
	  dati = datExpr[!del, i]; datj = datExpr[!del, j];
        lmPij=glm( dati ~ poly( datj, degree))
        polyRsquare[i, j] = cor( dati, predict(lmPij))^2
        lmPji=glm( datj ~ poly( dati, degree))
        polyRsquare[j, i] = cor( datj, predict(lmPji))^2
        rm(dati, datj, lmPij, lmPji)
      }
    }
  }

  diag(polyRsquare) = rep(1,n)

if (symmetrizationMethod =="none") {adj= polyRsquare} else
 {  adj = switch(symmetrizationMethod, 
	min = pmin(polyRsquare, t(polyRsquare)),
	max = pmax(polyRsquare, t(polyRsquare)),
	mean = (polyRsquare + t(polyRsquare))/2)
}
  adj

}


