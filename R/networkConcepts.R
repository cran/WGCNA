# ===================================================
# This code was written by Jun Dong, modified by Peter Langfelder
# datExpr: expression profiles with rows=samples and cols=genes/probesets
# power: for contruction of the weighted network
# trait: the quantitative external trait
#

networkConcepts = function(datExpr, power=1, trait=NULL, networkType = "unsigned") 
{

  networkTypeC = charmatch(networkType, .networkTypes);
  if (is.na(networkTypeC))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")));

  if(networkTypeC==1)
  {
        adj <- abs(cor(datExpr,use="p"))^power
  } else if (networkTypeC==2)
  {
  	adj <- abs((cor(datExpr,use="p")+1)/2)^power
  } else {
        cor = cor(datExpr,use="p");
        cor[cor < 0] = 0;
  	adj <- cor^power
  }
  diag(adj)=0 # Therefore adj=A-I.
	
  ### Fundamental Network Concepts
  Size=dim(adj)[1]
  Connectivity=apply(adj, 2, sum) # Within Module Connectivities
  Density=sum(Connectivity)/(Size*(Size-1))
  Centralization=Size*(max(Connectivity)-mean(Connectivity))/((Size-1)*(Size-2))
  Heterogeneity=sqrt(Size*sum(Connectivity^2)/sum(Connectivity)^2-1)
  ClusterCoef=.ClusterCoef.fun(adj)
  
  fMAR=function(v) sum(v^2)/sum(v)
  MAR=apply(adj, 1, fMAR)
  #CONNECTIVITY=Connectivity/max(Connectivity)
  
  ### Conformity-Based Network Concepts	
  ### Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24
  Conformity=.NPC.iterate(adj)$v1
  Factorizability=1- sum( (adj-outer(Conformity,Conformity)+ diag(Conformity^2))^2 )/sum(adj^2)
  Connectivity.CF=sum(Conformity)*Conformity-Conformity^2
  Density.CF=sum(Connectivity.CF)/(Size*(Size-1))
  Centralization.CF=Size*(max(Connectivity.CF)-mean(Connectivity.CF))/((Size-1)*(Size-2))
  Heterogeneity.CF=sqrt(Size*sum(Connectivity.CF^2)/sum(Connectivity.CF)^2-1)
  #ClusterCoef.CF=.ClusterCoef.fun(outer(Conformity,Conformity)-diag(Conformity^2) )
  ClusterCoef.CF=c(NA, Size)
  for(i in 1:Size )
     ClusterCoef.CF[i]=( sum(Conformity[-i]^2)^2 - sum(Conformity[-i]^4) )/
       ( sum(Conformity[-i])^2 - sum(Conformity[-i]^2) )

  ### Approximate Conformity-Based Network Concepts	
  Connectivity.CF.App=sum(Conformity)*Conformity
  Density.CF.App=sum(Connectivity.CF.App)/(Size*(Size-1))
  Centralization.CF.App=Size*(max(Connectivity.CF.App)-mean(Connectivity.CF.App))/((Size-1)*(Size-2))
  Heterogeneity.CF.App=sqrt(Size*sum(Connectivity.CF.App^2)/sum(Connectivity.CF.App)^2-1)
  ClusterCoef.CF.App=(sum(Conformity^2)/sum(Conformity))^2

  ### Eigengene-based Network Concepts
  m1=moduleEigengenes(datExpr, colors = rep(1, Size));
  # Weighted Expression Conformity
  ConformityE=cor(datExpr,m1[[1]][,1],use="pairwise.complete.obs"); ConformityE=abs(ConformityE)^power; 
  ConnectivityE=sum(ConformityE)*ConformityE; #Expression Connectivity
  DensityE=sum(ConnectivityE)/(Size*(Size-1)); #Expression Density
  CentralizationE=Size*(max(ConnectivityE)-mean(ConnectivityE))/((Size-1)*(Size-2)); #Expression Centralization
  HeterogeneityE=sqrt(Size*sum(ConnectivityE^2)/sum(ConnectivityE)^2-1); #Expression Heterogeneity
  ClusterCoefE=(sum(ConformityE^2)/sum(ConformityE))^2; ##Expression ClusterCoef
  MARE=ConformityE* sum(ConformityE^2)/sum(ConformityE)

  ### Significance measure only when trait is available.
  if(!is.null(trait)){
    EigengeneSignificance = abs(cor(trait, m1[[1]], use="pairwise.complete.obs") )^power;
    EigengeneSignificance = EigengeneSignificance[1,1]
    GS= abs(cor(datExpr, trait, use="pairwise.complete.obs") )^power; GS=GS[,1]
    GSE=ConformityE * EigengeneSignificance; GSE=GSE[,1]
    ModuleSignificance=mean(GS)
    ModuleSignificanceE=mean(GSE)
    K=Connectivity/max(Connectivity)
    HubGeneSignificance=sum(GS*K)/sum(K^2)
    KE=ConnectivityE/max(ConnectivityE)
    HubGeneSignificanceE= sum(GSE*KE)/sum(KE^2)
  }

  Summary=cbind(
    c(Density, Centralization, Heterogeneity, mean(ClusterCoef), mean(Connectivity)),
    c(DensityE, CentralizationE, HeterogeneityE, mean(ClusterCoefE), mean(ConnectivityE)),
    c(Density.CF, Centralization.CF, Heterogeneity.CF, mean(ClusterCoef.CF), mean(Connectivity.CF)),
    c(Density.CF.App, Centralization.CF.App, Heterogeneity.CF.App, mean(ClusterCoef.CF.App),
mean(Connectivity.CF.App) ) )
  colnames(Summary)=c("Fundamental", "Eigengene-based", "Conformity-Based", "Approximate Conformity-based")
  rownames(Summary)=c("Density", "Centralization", "Heterogeneity", "Mean ClusterCoef", "Mean Connectivity")

  output=list(Summary=Summary, Size=Size, Factorizability=Factorizability, Eigengene=m1[[1]],
VarExplained=m1[[2]][,1], Conformity=Conformity, ClusterCoef=ClusterCoef, Connectivity=Connectivity,
MAR=MAR, ConformityE=ConformityE)
  if(!is.null(trait)){
    output$GS=GS; output$GSE=GSE;
    Significance=cbind(c(ModuleSignificance, HubGeneSignificance, EigengeneSignificance),
    c(ModuleSignificanceE, HubGeneSignificanceE, NA))
    colnames(Significance)=c("Fundamental", "Eigengene-based")
    rownames(Significance)=c("ModuleSignificance", "HubGeneSignificance", "EigengeneSignificance")
    output$Significance=Significance
  }
	output
}



#====================================================================================================
#
# Network functions for network concepts
#
#====================================================================================================
#=========================================
# Function definitions
#=========================================

# ================================================================================
#   Cohesiveness/Conformity/Factorizability etc
# ================================================================================

# ===================================================
# Check if adj is a valid adjacency matrix:  square matrix, non-negative entries, symmetric and no missing entries.
# Parameters: 
#   adj - the input adjacency matrix
#   tol - the tolerence level to measure the difference from 0 (symmetric matrix: upper diagonal minus lower diagonal)
# Remarks:
#   1. This function is not supposed to be used directly. Instead, it should appear in function definitions.
#   2. We release the requirement that the diagonal elements be 1 or 0. Users should assign appropriate values
#      at the beginning of their function definitions.
# Usage:
# if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")

.is.adjmat = function(adj, tol=10^(-15)){
	n=dim(adj)
	is.adj=1
	if (n[1] != n[2]){ message("The adjacency matrix is not a square matrix!"); is.adj=0;}
	if ( sum(is.na(adj))>0 ){ message("There are missing values in the adjacency matrix!"); is.adj=0;}
	if ( sum(adj<0)>0 ){ message("There are negative entries in the adjacency matrix!"); is.adj=0;}
	if ( max(abs(adj-t(adj))) > tol){ message("The adjacency matrix is not symmetric!"); is.adj=0;}
			#if ( max(abs(diag(adj)-1)) > tol){ message("The diagonal elements are not all one!"); is.adj=0;}
			#The last criteria is removed because of different definitions on diagonals with other papers.
			#Always let "diagonal=1" INSIDE the function calls when using functions for Factorizability paper.
	is.adj
}


# ===================================================
# .NPC.direct=function(adj)
# Calculates the square root of Normalized Product Connectivity (.NPC), by way of definition of .NPC. ( \sqrt{t})
# Parameters: 
#   adj - the input adjacency matrix
#   tol - the tolerence level to measure the difference from 0 (zero off-diagonal elements)
# Output:
#   v1 - vector, the square root of .NPC
# Remarks:
#   1. The function requires that the off-diagonal elements of the adjacency matrix are all non-zero.
#   2. If any of the off-diagonal elements is zero, use the function .NPC.iterate().
#   3. If the adjacency matrix is 2 by 2, then a warning message is issued and vector of sqrt(adj[1,2]) is returned.
#   4. If the adjacency matrix is a ZERO matrix, then a warning message is issued and vector of 0 is returned.

.NPC.direct=function(adj){
	if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
	n=dim(adj)[1]
	if(n==2) {
		warning("The adjacecny matrix is only 2 by 2. .NPC may not be unique!")
		return(rep(sqrt(adj[1,2]),2))
	}
	diag(adj)=0
	if(!sum(adj>0)){
	  warning("The adjacency matrix is a ZERO matrix!")
	  return(rep(0,n))
	}
	diag(adj)=1
	if(sum(adj==0)) stop("There is zero off--diagonal element! Please use the function .NPC.iterate().")
	log10.prod.vec=function(vec){
		prod=0
		for(i in 1:length(vec) )
			prod=prod+log10(vec[i])
		prod
	}
	off.diag=as.vector(as.dist(adj))
	prod1=log10.prod.vec(off.diag)
	v1=rep(-666, n)
	for(i in 1:n){
		prod2=prod1-log10.prod.vec(adj[i,])
		v1[i]=10^(prod1/(n-1)-prod2/(n-2))
	}
	v1
}

# ===================================================
# .NPC.iterate=function(adj, loop=10^(10), tol=10^(-10))
# Calculates the square root of Normalized Product Connectivity, by way of iteration algorithm. ( \sqrt{t})
# Parameters: 
#   adj - the input adjacency matrix
#   loop - the maximum number of iterations before stopping the algorithm
#   tol - the tolerence level to measure the difference from 0 (zero off-diagonal elements)
# Output:
#   v1 - vector, the square root of .NPC
#   loop - integer, the number of iterations taken before convergence criterion is met
#   diff - scaler, the maximum difference between the estimates of 'v1' in the last two iterations
# Remarks:
#   1. Whenever possible, use .NPC.direct().
#   2. If the adjacency matrix is 2 by 2, then a warning message is issued.
#   3. If the adjacency matrix is a ZERO matrix, then a warning message is issued and vector of 0 is returned.

if( exists(".NPC.iterate") ) rm(.NPC.iterate);
.NPC.iterate=function(adj, loop=10^(10), tol=10^(-10)){
	if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
	n=dim(adj)[1]
	if(n==2) warning("The adjacecny matrix is only 2 by 2. .NPC may not be unique!")
	diag(adj)=0
	if(max(abs(adj))<tol){
	  warning("The adjacency matrix is a ZERO matrix!")
	  return(rep(0,n))
	}
	diff=1
	k=apply(adj, 2, sum)
	v1=k/sqrt(sum(k)) # first-step estimator for v-vector (initial value)
	i=0
	while( loop>i && diff>tol ){
	i=i+1
	diag(adj)=v1^2
	svd1=svd(adj) # Spectral Decomposition
	v2=sqrt(svd1$d[1])*abs(svd1$u[,1])
	diff=max(abs(v1-v2))
	v1=v2
	}
	list(v1=v1,loop=i,diff=diff)
}

# ===================================================
# The function .ClusterCoef.fun computes the cluster coefficients.
# Input is an adjacency matrix 
.ClusterCoef.fun=function(adjmat1) 
{
  # diag(adjmat1)=0
  no.nodes=dim(adjmat1)[[1]]
  computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
  computeSqDiagSum = function(x, vec) { sum(x^2 * vec) };
  nolinksNeighbors <- c(rep(-666,no.nodes))
  total.edge <- c(rep(-666,no.nodes))
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) 
  {
     stop(paste("ERROR: the adjacency matrix contains entries that are larger",
                 "than 1 or smaller than 0: max=",maxh1,", min=",minh1)) 
  } else { 
    nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, imatrix=adjmat1)
    subTerm = apply(adjmat1, 1, computeSqDiagSum, vec = diag(adjmat1));
    plainsum  <- apply(adjmat1, 1, sum)
    squaresum <- apply(adjmat1^2, 1, sum)
    total.edge = plainsum^2 - squaresum
    CChelp=rep(-666, no.nodes)
    CChelp=ifelse(total.edge==0,0, (nolinksNeighbors-subTerm)/total.edge)
    CChelp
  }
} # end of function


# ===================================================
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)
.err.bp<-function(daten,error,two.side=F)
{
   if(!is.numeric(daten)) {
        stop("All arguments must be numeric")}
   if(is.vector(daten)){ 
      xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
   }else{
      if (is.matrix(daten)){
        xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
                            dim=c(1,length(daten))))+0:(length(daten)-1)+.5
      }else{
        stop("First argument must either be a vector or a matrix") }
   }
   MW<-0.25*(max(xval)/length(xval)) 
   ERR1<-daten+error 
   ERR2<-daten-error
   for(i in 1:length(daten)){
      segments(xval[i],daten[i],xval[i],ERR1[i])
      segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
      if(two.side){
        segments(xval[i],daten[i],xval[i],ERR2[i])
        segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
      } 
   } 
} 



#========================================================================================

conformityBasedNetworkConcepts = function(adj, GS=NULL) 
{
  if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
  diag(adj)=0 # Therefore adj=A-I.	
  if (dim(adj)[[1]]<3) 
    stop("The adjacency matrix has fewer than 3 rows. This network is trivial and will not be evaluated.")
  if (!is.null(GS)) 
  { 
    if( length(GS) !=dim(adj)[[1]])
    { 
       stop(paste("The length of the node significnce GS does not equal the number",
                  "of rows of the adjcency matrix. length(GS) != dim(adj)[[1]]. \n",
                  "Something is wrong with your input"))
    }
  }

	### Fundamental Network Concepts
	Size=dim(adj)[1]
	Connectivity=apply(adj, 2, sum) 
	Density=sum(Connectivity)/(Size*(Size-1))
	Centralization=Size*(max(Connectivity)-mean(Connectivity))/((Size-1)*(Size-2))
	Heterogeneity=sqrt(Size*sum(Connectivity^2)/sum(Connectivity)^2-1)
	ClusterCoef=.ClusterCoef.fun(adj)
	fMAR=function(v) sum(v^2)/sum(v)
	MAR=apply(adj, 1, fMAR)	
	### Conformity-Based Network Concepts	
	Conformity=.NPC.iterate(adj)$v1
	Factorizability=1- sum( (adj-outer(Conformity,Conformity)+ diag(Conformity^2))^2 )/sum(adj^2)
	Connectivity.CF=sum(Conformity)*Conformity-Conformity^2
	Density.CF=sum(Connectivity.CF)/(Size*(Size-1))
	Centralization.CF=Size*(max(Connectivity.CF)-mean(Connectivity.CF))/((Size-1)*(Size-2))
	Heterogeneity.CF=sqrt(Size*sum(Connectivity.CF^2)/sum(Connectivity.CF)^2-1)
	#ClusterCoef.CF=.ClusterCoef.fun(outer(Conformity,Conformity)-diag(Conformity^2) )
	ClusterCoef.CF=c(NA, Size)
  for(i in 1:Size )
    ClusterCoef.CF[i]=( sum(Conformity[-i]^2)^2 - sum(Conformity[-i]^4) )/
                          ( sum(Conformity[-i])^2 - sum(Conformity[-i]^2) )
  MAR.CF=ifelse(sum(Conformity,na.rm=T)-Conformity==0, NA,
                Conformity*(sum(Conformity^2,na.rm=T)-Conformity^2)/(sum(Conformity,na.rm=T)-Conformity))

	### Approximate Conformity-Based Network Concepts	
	Connectivity.CF.App=sum(Conformity)*Conformity
	Density.CF.App=sum(Connectivity.CF.App)/(Size*(Size-1))
	Centralization.CF.App=Size*(max(Connectivity.CF.App)-mean(Connectivity.CF.App))/((Size-1)*(Size-2))
	Heterogeneity.CF.App=sqrt(Size*sum(Connectivity.CF.App^2)/sum(Connectivity.CF.App)^2-1)
	
  if(sum(Conformity,na.rm=T)==0)
  {
      warning(paste("The sum of conformities equals zero.\n",
                    "Maybe you used an input adjacency matrix with lots of zeroes?\n",
                    "Specifically, sum(Conformity,na.rm=T)==0."));
      MAR.CF.App= rep(NA,Size) 
      ClusterCoef.CF.App= rep(NA,Size) 
  } #end of if
  if(sum(Conformity,na.rm=T) !=0)
  {
    MAR.CF.App=Conformity*sum(Conformity^2,na.rm=T) /sum(Conformity,na.rm=T)
    ClusterCoef.CF.App=rep((sum(Conformity^2)/sum(Conformity))^2,Size)
  }# end of if
  output=list(
              Factorizability =Factorizability, 
              fundamentalNCs=list(
                  ScaledConnectivity=Connectivity/max(Connectivity,na.rm=T),
                  Connectivity=Connectivity, 
                  ClusterCoef=ClusterCoef, 
                  MAR=MAR, 
                  Density=Density,
                  Centralization =Centralization,
                  Heterogeneity= Heterogeneity),
              conformityBasedNCs=list(
                  Conformity=Conformity, 
                  Connectivity.CF=Connectivity.CF, 
                  ClusterCoef.CF=ClusterCoef.CF,
                  MAR.CF=MAR.CF,
                  Density.CF=Density.CF,
                  Centralization.CF =Centralization.CF,
                  Heterogeneity.CF= Heterogeneity.CF),
              approximateConformityBasedNCs=list(
                  Conformity=Conformity, 
                  Connectivity.CF.App= Connectivity.CF.App,
                  ClusterCoef.CF.App=ClusterCoef.CF.App,
                  MAR.CF.App=MAR.CF.App,
                  Density.CF.App= Density.CF.App,
                  Centralization.CF.App =Centralization.CF.App,
                  Heterogeneity.CF.App= Heterogeneity.CF.App))
  if ( !is.null(GS) ) 
  {
    output$FundamentalNC$NetworkSignificance = mean(GS,na.rm=T)
    K = Connectivity/max(Connectivity)
    output$FundamentalNC$HubNodeSignificance = sum(GS * K,na.rm=T)/sum(K^2,na.rm=T)	
  }
  output
} # end of function 



#===================================================================================================

fundamentalNetworkConcepts=function(adj,GS=NULL)
{
   if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
   diag(adj)=0 # Therefore adj=A-I.
	
   if (dim(adj)[[1]]<3) stop("The adjacency matrix has fewer than 3 rows. This network is trivial and will not be evaluated.")

if (!is.null(GS)) { if( length(GS) !=dim(adj)[[1]]){ stop("The length of the node significnce GS does not
equal the number of rows of the adjcency matrix. length(GS) unequal dim(adj)[[1]]. GS should be a vector whosecomponents correspond to the nodes.")}}

	Size=dim(adj)[1]
	
### Fundamental Network Concepts
	Connectivity=apply(adj, 2, sum) # Within Module Connectivities
	Density=sum(Connectivity)/(Size*(Size-1))
	Centralization=Size*(max(Connectivity)-mean(Connectivity))/((Size-1)*(Size-2))
	Heterogeneity=sqrt(Size*sum(Connectivity^2)/sum(Connectivity)^2-1)
	ClusterCoef=.ClusterCoef.fun(adj)
	fMAR=function(v) sum(v^2)/sum(v)
	MAR=apply(adj, 1, fMAR)
	ScaledConnectivity=Connectivity/max(Connectivity,na.rm=T)
	
	output=list(
                  Connectivity=Connectivity, 
                  ScaledConnectivity=ScaledConnectivity, 
                  ClusterCoef=ClusterCoef, MAR=MAR, 
                  Density=Density, Centralization =Centralization,
                  Heterogeneity= Heterogeneity) 

   if ( !is.null(GS) ) {
        output$NetworkSignificance = mean(GS,na.rm=T)
        output$HubNodeSignificance = sum(GS * ScaledConnectivity,na.rm=T)/sum(ScaledConnectivity^2,na.rm=T)
   }
   output
} # end of function 

#==========================================================================================================
#
# Density function
#
#==========================================================================================================




