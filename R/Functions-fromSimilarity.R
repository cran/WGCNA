# Functions to perform WGCNA from similarity input.

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
  pickHardThreshold(similarity,  RsquaredCut =  RsquaredCut, cutVector = cutVector,
                    moreNetworkConcepts = moreNetworkConcepts, removeFirst = removeFirst,
                    nBreaks = nBreaks, corFnc = "I", corOptions = "");
}

pickSoftThreshold.fromSimilarity = function (similarity, 
    RsquaredCut = 0.85, powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
    removeFirst = FALSE, nBreaks = 10, blockSize = 1000,
    networkType = "unsigned", moreNetworkConcepts=FALSE, verbose = 0, indent = 0)
{
  checkSimilarity(similarity)
  pickSoftThreshold(similarity,  RsquaredCut =  RsquaredCut, powerVector = powerVector,
                    removeFirst = removeFirst, nBreaks = nBreaks, 
                    blockSize = blockSize, networkType = networkType,
                    moreNetworkConcepts = moreNetworkConcepts,
                    verbose = verbose, indent = indent);

}


