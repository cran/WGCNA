# Definitions of internal constants.

.pearsonFallbacks = c("none", "individual", "all");

.zeroMADWarnings = c("Some results will be NA.", 
                     "Pearson correlation was used for individual columns with MAD=NA.",
                     "Pearson correlation was used for entire variable.");

..minNGenes = 4;
..minNSamples = 4;

.largestBlockSize = 1e8;

.networkTypes = c("unsigned", "signed", "signed hybrid");
.adjacencyTypes = c(.networkTypes, "distance");
.TOMTypes = c("none", "unsigned", "signed");

.TOMDenoms = c("min", "mean");

.corTypes = c("pearson", "bicor");

.corFnc = c("cor", "bicor", "cor");
.corOptions = c("use = 'p'", "use = 'p'", "use = 'p', method = spearman");


