# Functions to perform WGCNA from similarity input.

matrixToNetwork = function(mat,
   symmetrizeMethod = c("average", "min", "max"),
   signed = TRUE,
   min = NULL,
   max = NULL,
   power = 12,
   diagEntry = 1)
{
  sm = match.arg(symmetrizeMethod);
  if (is.na(sm)) 
    stop("Unrecognized or non-unique 'symmetrizeMethod'.");

  mat = as.matrix(mat);

  nd = 0
  x = try({nd = dim(mat)});
  if ( (class(x)=='try-error') | (nd!=2) )
    stop("'mat' appears to have incorrect type; must be a 2-dimensional square matrix.");

  if (ncol(mat)!=nrow(mat))
    stop("'mat' must be a square matrix.");

  if (!signed) mat = abs(mat);

  if (sm==1) {
    mat = (mat + t(mat))/2;
  } else if (sm==2) {
    mat = pmin(mat, t(mat), na.rm = TRUE);
  } else
    mat = pmax(mat, t(mat), na.rm = TRUE);
 
  if (is.null(min)) {
    min = min(mat, na.rm = TRUE);
  } else
    mat[mat < min] = min;

  if (is.null(max)) {
    max = max(mat, na.rm = TRUE);
  } else 
    mat[mat > max] = max;

  adj = ( (mat-min)/(max-min) )^power;

  diag(adj) = diagEntry

  adj;
}

checkSimilarity = function(similarity, min=-1, max=1)
{
  checkAdjMat(similarity, min, max);
}

adjacency.fromSimilarity = function(similarity, type = "unsigned", power = if (type=="distance") 1 else 6)
{
  checkSimilarity(similarity);
  adjacency(similarity, type = type, power = power, corFnc = "I", corOptions="", distFnc = "I", 
            distOptions = "");
}

softConnectivity.fromSimilarity=function(similarity, type = "unsigned",
                          power = if (type == "signed") 15 else 6,
                          blockSize = 1500, verbose = 2, indent = 0)

{
  checkSimilarity(similarity)
  softConnectivity(similarity, corFnc = "I", corOptions = "", 
                   type = type, power = power, 
                   blockSize = blockSize, verbose = verbose, indent = indent)
}

pickHardThreshold.fromSimilarity=function (similarity,
    RsquaredCut = 0.85, cutVector = seq(0.1, 0.9, by = 0.05), 
    moreNetworkConcepts=FALSE , removeFirst = FALSE, nBreaks = 10)
{
  checkSimilarity(similarity)
  pickHardThreshold(similarity, dataIsExpr = FALSE, RsquaredCut =  RsquaredCut, cutVector = cutVector,
                    moreNetworkConcepts = moreNetworkConcepts, removeFirst = removeFirst,
                    nBreaks = nBreaks, corFnc = "I", corOptions = "");
}

pickSoftThreshold.fromSimilarity = function (similarity, 
    RsquaredCut = 0.85, powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
    removeFirst = FALSE, nBreaks = 10, blockSize = 1000,
    networkType = "unsigned", moreNetworkConcepts=FALSE, verbose = 0, indent = 0)
{
  checkSimilarity(similarity)
  pickSoftThreshold(similarity, dataIsExpr = FALSE, RsquaredCut =  RsquaredCut, powerVector = powerVector,
                    removeFirst = removeFirst, nBreaks = nBreaks, 
                    blockSize = blockSize, networkType = networkType,
                    moreNetworkConcepts = moreNetworkConcepts,
                    verbose = verbose, indent = indent);

}


