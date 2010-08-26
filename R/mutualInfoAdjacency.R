#read in mutual information fx


mutualInfoAdjacency=function(  datE, discretizeColumns=TRUE, entropyEstimationMethod="MM", numberBins=NULL) {

reqPackages = c("infotheo", "minet", "entropy")
nReq = length(reqPackages)

for (r in 1:nReq)
{
  expression = spaste('require("', reqPackages[r], '", quietly = TRUE)');
  ok = eval(parse(text = expression))
  if (!ok) stop("The function requires R package infotheo, minet and entropy. Please install these packages first.")
}

if ( !is.element( discretizeColumns, c(TRUE, FALSE) ) ) stop("The input parameter discretizeColumns
contains a value that is not logical. It needs to be set to TRUE or FALSE")
if ( !is.element( entropyEstimationMethod, c("MM", "ML", "shrink", "SG" ) ) ){warning("The entropy estimation method does not correspond to any of the following: MM, ML, shrink, SG.  MM will be used."); entropyEstimationMethod="MM"}
datE=data.frame(datE)
if ( ! (dim(datE)[[2]]>1) ) stop("The number of columns of datE must be larger than 1")

entropyOfCountData=function( counts ) {
  express='entropy::entropy(table(counts) , unit="log", method= entropyEstimationMethod)'
  eval(parse(text=express))
}
if (is.null(numberBins) )   numberBins=sqrt(nrow(datE))   
if (!is.null(numberBins) )   numberBins=as.integer(numberBins)
 if ( !( numberBins>1 ) ) stop("Something is wrong with the input parameter numberBins, which is used for discretizing the quantitative variables. numberBins should be larger than 1. Recommendation: choose the default value numberBins=NULL")  

if (discretizeColumns) { 
  express = 'infotheo::discretize(datE, disc ="equalwidth", nbins=numberBins)'
  discretized.datE= eval(parse(text = express));
}
if (! discretizeColumns) { discretized.datE= datE }
ENTROPY=as.numeric( apply(discretized.datE, 2, entropyOfCountData ))
if (entropyEstimationMethod =="MM") entropyEstimationMethodRenamed="mi.mm" ;
if (entropyEstimationMethod =="ML") entropyEstimationMethodRenamed="mi.empirical" ;
if (entropyEstimationMethod =="shrink") entropyEstimationMethodRenamed="mi.shrink" ;
if (entropyEstimationMethod =="SG") entropyEstimationMethodRenamed="mi.sg" ;

express="minet::build.mim(discretized.datE , estimator= entropyEstimationMethodRenamed)"
MIxy = eval(parse(text=express))
MIxy[MIxy<0]=0
diag(MIxy)=ENTROPY
AdjacencySymmetricUncertainty=2*MIxy/ outer(ENTROPY,ENTROPY, FUN="+")  
AdjacencyUniversal= AdjacencySymmetricUncertainty/(2- AdjacencySymmetricUncertainty)
AdjacencyUniversalVersion2= MIxy/outer(ENTROPY,ENTROPY, FUN="pmax", na.rm=T)  
#output
list(Entropy=ENTROPY, MutualInformation=MIxy, 
AdjacencySymmetricUncertainty= AdjacencySymmetricUncertainty, AdjacencyUniversalVersion1= AdjacencyUniversal, AdjacencyUniversalVersion2= AdjacencyUniversalVersion2)
} # end of function mutualInfoAdjacency

