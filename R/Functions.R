# Categories of functions:
# . network construction (including connectivity calculation)
# . module detection
# . gene screening
# . data simulation
# . general statistical functions
# . visualization


#==========================================================================================
#
#  Subselection for multi-set data structure
#
#==========================================================================================

#subselect(multiData, index, selectOnComps = 1)
#{
#  size = checkSets(multiData);
#  selected = multiData;
#  for (set in 1:size$nSets)
#  {
#    for (comp in selectOnComps)
#    {
#      selected[[set]][[comp]] = multiData[[set]][[comp]][,index]
#    }
#  }
#}
  

# ===================================================
#For hard thresholding, we use the signum (step) function

signumAdjacencyFunction = function(corMat, threshold)  
{
  adjmat= as.matrix(abs(corMat)>=threshold)
  dimnames(adjmat) <- dimnames(corMat)
  diag(adjmat) <- 0
  adjmat
}

# ===================================================
# For soft thresholding, one can use the sigmoid function 
# But we have focused on the power adjacency function in the tutorial...
sigmoidAdjacencyFunction = function(ss, mu=0.8, alpha=20) 
{
   1/(1+exp(-alpha*(ss-mu)))
}

#This function is useful for speeding up the connectivity calculation.
#The idea is to partition the adjacency matrix into consecutive baches of a #given size.
#In principle, the larger the block size the faster is the calculation. But #smaller blockSizes require #less memory...
# Input: gene expression data set where *rows* correspond to microarray samples #and columns correspond to genes. 
# If fewer than minNSamples contain gene expression information for a given
# gene, then its connectivity is set to missing. 

softConnectivity=function(datExpr, 
                          corFnc = "cor", corOptions = "use = 'p'", 
                          type = "unsigned", 
                          power = if (type == "signed") 15 else 6, 
                          blockSize = 1500, minNSamples = NULL,
                          verbose = 2, indent = 0) 
{
  spaces = indentSpaces(indent);
  nGenes=dim(datExpr)[[2]]

  if (blockSize * nGenes>.largestBlockSize) blockSize = as.integer(.largestBlockSize/nGenes);
  nSamples=dim(datExpr)[[1]]
  if (is.null(minNSamples))
  {
    minNSamples = max(..minNSamples, nSamples/3);
  }

  if (nGenes<..minNGenes | nSamples<minNSamples ) 
     stop(paste("Error: Something seems to be wrong. \n", 
          "   Make sure that the input data frame has genes as rows and array samples as columns.\n",
          "   Alternatively, there seem to be fewer than", ..minNGenes, "genes or fewer than", 
              minNSamples, "samples. ") )
  if (nGenes<nSamples ) 
    printFlush("Warning: There are fewer genes than samples in the function softConnectivity. Maybe you should transpose the data?")


  k=rep(NA,nGenes)
  start = 1;
  if (verbose>0) { 
     printFlush(paste(spaces, "softConnectivity: FYI: connecitivty of genes with less than", 
                               ceiling(minNSamples), "valid samples will be returned as NA.")); 
     cat(paste(spaces, "..calculating connectivities..")); 
     pind = initProgInd();
  }
  while (start < nGenes)
  {
    end = min(start + blockSize-1, nGenes);
    index1=start:end;
    ad1 = adjacency(datExpr, index1, power = power, type = type, 
                    corFnc = corFnc, corOptions = corOptions);
    k[index1]=colSums(ad1, na.rm = TRUE)-1;
    # If fewer than minNSamples contain gene expression information for a given
    # gene, then we set its connectivity to 0.
    NoSamplesAvailable=colSums(!is.na(datExpr[,index1]))
    k[index1][NoSamplesAvailable< minNSamples]=NA
    if (verbose>0) pind = updateProgInd(end/nGenes, pind);
    start = end + 1;
  } 
  if (verbose > 0) printFlush("");
  k
} # end of function


# ==============================================================================
# The function PickHardThreshold can help one to estimate the cut-off value 
# when using the signum (step) function.
# The first column lists the threshold ("cut"),
# the second column lists the corresponding p-value based on the Fisher Transform 
# of the correlation. 
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative for 
# biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model. 
# Usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative slope 
# (around -1, col 4).
# The output is a list with 2 components. The first component lists a sugggested cut-off
# while the second component contains the whole table.
# The removeFirst option removes the first point (k=0, P(k=0)) from the regression fit.
# nBreaks specifies how many intervals used to estimate the frequency p(k) i.e. the no. of points in the 
# scale free topology plot.

pickHardThreshold=function (data, dataIsExpr = TRUE, RsquaredCut = 0.85, cutVector = seq(0.1, 0.9, 
    by = 0.05), moreNetworkConcepts=FALSE , removeFirst = FALSE, nBreaks = 10, corFnc = "cor", 
    corOptions = "use = 'p'") 
{
    nGenes = dim(data)[[2]]
    colname1 = c("Cut", "p-value", "SFT.R.sq", "slope=", 
        "truncated R^2", "mean(k)", "median(k)", "max(k)")
if(moreNetworkConcepts) {
colname1=c(colname1,"Density", "Centralization", "Heterogeneity")
}
    if (!dataIsExpr)
    {
      checkAdjMat(data);
      if (any(diag(data)!=1)) diag(data) = 1;
    } else
      nSamples = dim(data)[[1]]

    datout = data.frame(matrix(NA, nrow = length(cutVector), 
        ncol = length(colname1)))
    names(datout) = colname1
    datout[, 1] = cutVector
    for (i in 1:length(cutVector)) {
        cut1 = cutVector[i]
        datout[i, 2] = 2 * (1 - pt(sqrt(nSamples - 1) * cut1/sqrt(1 - 
            cut1^2), nSamples - 1))
    }
    if (exists("fun1")) rm(fun1)
    fun1 = function(x, dataIsExpr) {
        if (dataIsExpr)
        {
          corExpr = parse(text = paste(corFnc, "(x, data", 
              prepComma(corOptions), ")"))
          corx = abs(eval(corExpr))
        } else 
          corx = x;
        out1 = rep(NA, length(cutVector))
        for (j in c(1:length(cutVector))) {
            out1[j] = sum(corx > cutVector[j])
        }
        out1
    }
    datk = t(apply(data, 2, fun1, dataIsExpr))
    for (i in c(1:length(cutVector))) {
        khelp= datk[, i] - 1
        SFT1=scaleFreeFitIndex(k=khelp,nBreaks=nBreaks,removeFirst=removeFirst)
        datout[i, 3] = SFT1$Rsquared.SFT  
        datout[i, 4] = SFT1$slope.SFT 
        datout[i, 5] = SFT1$truncatedExponentialAdjRsquared
        datout[i, 6] = mean(khelp,na.rm=T)
        datout[i, 7] = median(khelp,na.rm=T)
        datout[i, 8] = max(khelp,na.rm=T)
if(moreNetworkConcepts) { 
Density = sum(khelp)/(nGenes * (nGenes - 1))
datout[i, 9] =Density
Centralization = nGenes*(max(khelp)-mean(khelp))/((nGenes-1)*(nGenes-2))
datout[i, 10] = Centralization
Heterogeneity = sqrt(nGenes * sum(khelp^2)/sum(khelp)^2 - 1)
datout[i, 11] = Heterogeneity
}
    }
    print(signif(data.frame(datout),3))
    ind1 = datout[, 3] > RsquaredCut
    indcut = NA
    indcut = ifelse(sum(ind1) > 0, min(c(1:length(ind1))[ind1]), 
        indcut)
    cutEstimate = cutVector[indcut][[1]]
    list(cutEstimate = cutEstimate, fitIndices = data.frame(datout))
} # end of function pickHardThreshold


#==============================================================================================
#
# pickSoftThreshold
#
#===============================================================================================
# The function pickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
# PL: a rewrite that splits the data into a few blocks.
# SH: more netowkr concepts added.
# PL: re-written for parallel processing

pickSoftThreshold = function (data, dataIsExpr = TRUE, RsquaredCut = 0.85, 
    powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), removeFirst = FALSE, nBreaks = 10, 
    blockSize = NULL, corFnc = cor, corOptions = list(use = 'p'), 
    networkType = "unsigned", moreNetworkConcepts=FALSE, verbose = 0, indent = 0)
{
    intType = charmatch(networkType, .networkTypes)
    if (is.na(intType)) 
        stop(paste("Unrecognized 'networkType'. Recognized values are", 
            paste(.networkTypes, collapse = ", ")))
    nGenes = ncol(data);
    if (nGenes<3) 
    { 
       stop("The input data data contain fewer than 3 rows (nodes).", 
            "\nThis would result in a trivial correlation network." )
    }
    if (!dataIsExpr) 
    {
      checkAdjMat(data);
      if (any(diag(data)!=1)) diag(data) = 1;
    }

    if (is.null(blockSize))
    {
      blockSize = blockSize(nGenes, rectangularBlocks = TRUE);
      if (verbose > 0) 
        printFlush(spaste("pickSoftThreshold: will use block size ", blockSize, "."))
    }
   
    colname1 = c("Power", "SFT.R.sq", "slope", "truncated R.sq", 
                 "mean(k)", "median(k)", "max(k)")
    if(moreNetworkConcepts) 
    {
         colname1=c(colname1,"Density", "Centralization", "Heterogeneity")
    }
    datout = data.frame(matrix(666, nrow = length(powerVector), ncol = length(colname1)))
    names(datout) = colname1
    datout[, 1] = powerVector
    spaces = indentSpaces(indent)
    if (verbose > 0) {
        cat(paste(spaces, "pickSoftThreshold: calculating connectivity for given powers..."))
        if (verbose == 1) 
            pind = initProgInd()
        else cat("\n")
    }

    # if we're using one of WGNCA's own correlation functions, set the number of threads to 1.
    corFnc = match.fun(corFnc);
    corFormals = formals(corFnc);
    if ("nThreads" %in% names(corFormals)) corOptions$nThreads = 1;

    # Resulting connectivities
    datk = matrix(0, nrow = nGenes, ncol = length(powerVector))

    # Number of threads. In this case I need this explicitly.
    nThreads = WGCNAnThreads();

    nPowers = length(powerVector);

    # Main loop
    startG = 1
    while (startG <= nGenes) 
    {
      endG = min (startG + blockSize - 1, nGenes)

      if (verbose > 1) 
          printFlush(paste(spaces, "  ..working on genes", 
              startG, "through", endG, "of ", nGenes))

      nBlockGenes = endG - startG + 1;
      jobs = allocateJobs(nBlockGenes, nThreads);
      # This assumes that the non-zero length allocations
      # precede the zero-length ones
      actualThreads = which(sapply(jobs, length) > 0); 

      datk[ c(startG:endG), ] = foreach(t = actualThreads, .combine = rbind) %dopar% 
      {
        useGenes = c(startG:endG)[ jobs[[t]] ]
        nGenes1 = length(useGenes);
        if (dataIsExpr)
        {
          corOptions$x = data;
          corOptions$y = data[ , useGenes];
          corx = do.call(corFnc, corOptions);
          if (intType == 1) {
              corx = abs(corx)
          } else if (intType == 2) {
              corx = (1 + corx)/2
          } else if (intType == 3) {
              corx[corx < 0] = 0
          }
          if (sum(is.na(corx)) != 0) 
              warning(paste("Some correlations are NA in block", 
                  startG, ":", endG, "."))
        } else {
          corx = data[, useGenes];
        }
        datk.local = matrix(NA, nGenes1, nPowers);
        for (j in 1:nPowers) {
            datk.local[, j] = colSums(corx^powerVector[j], na.rm = TRUE) - 1;
        }
        datk.local;
      } # End of %dopar% evaluation
      # Move to the next block of genes.
      startG = endG + 1
      if (verbose == 1) pind = updateProgInd(endG/nGenes, pind)
    }
    if (verbose == 1) printFlush("");

    for (i in c(1:length(powerVector))) 
    {
        khelp= datk[, i] 
        SFT1=scaleFreeFitIndex(k=khelp,nBreaks=nBreaks,removeFirst=removeFirst)
        datout[i, 2] = SFT1$Rsquared.SFT  
        datout[i, 3] = SFT1$slope.SFT 
        datout[i, 4] = SFT1$truncatedExponentialAdjRsquared
        datout[i, 5] = mean(khelp,na.rm=T)
        datout[i, 6] = median(khelp,na.rm=T)
        datout[i, 7] = max(khelp,na.rm=T)
        if(moreNetworkConcepts) 
        { 
           Density = sum(khelp)/(nGenes * (nGenes - 1))
           datout[i, 8] =Density
           Centralization = nGenes*(max(khelp)-mean(khelp))/((nGenes-1)*(nGenes-2))
           datout[i, 9] = Centralization
           Heterogeneity = sqrt(nGenes * sum(khelp^2)/sum(khelp)^2 - 1)
           datout[i, 10] = Heterogeneity
         }
    }
    print(signif(data.frame(datout),3))
    ind1 = datout[, 2] > RsquaredCut
    indcut = NA
    indcut = ifelse(sum(ind1) > 0, min(c(1:length(ind1))[ind1]), indcut)
    powerEstimate = powerVector[indcut][[1]]
    list(powerEstimate = powerEstimate, fitIndices = data.frame(datout))
}


# ===================================================
# The function ScaleFreePlot1 creates a plot for checking scale free topology
# when truncated1=T is specificed, it provides the R^2 measures for the following
# degree distributions: a) scale free topology, b) log-log R^2 and c) truncated exponential R^2

# The function ScaleFreePlot1 creates a plot for checking scale free topology

scaleFreePlot = function(connectivity, nBreaks=10, truncated = FALSE, removeFirst = FALSE, main = "", ...)
{
  k = connectivity
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k),
      length = nBreaks + 1)
  hist1 = suppressWarnings(hist(k, breaks = breaks1, equidist = FALSE, plot = FALSE, right = TRUE, ...))
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
      p.dk = p.dk[-1]
      log.dk = log.dk[-1]
  }
  log.p.dk= as.numeric(log10(p.dk + 1e-09))
  lm1 = lm(log.p.dk ~ log.dk)
  if (truncated==TRUE) 
  { 
    lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
    OUTPUT=data.frame(scaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),
                      slope=round(lm1$coefficients[[2]],2),
    TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
    printFlush("the red line corresponds to the truncated exponential fit")
    title = paste(main, 
                " scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
                ", slope=", round(lm1$coefficients[[2]],2),
                ", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2)))
  } else { 
    title = paste(main, " scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
                  ", slope=", round(lm1$coefficients[[2]],2))
    OUTPUT=data.frame(scaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),
                      slope=round(lm1$coefficients[[2]],2))
  }

  suppressWarnings(plot(log.dk, log.p.dk, xlab="log10(k)", ylab="log10(p(k))", main = title, ... ))
  lines(log.dk,predict(lm1),col=1)
  if (truncated) lines(log.dk, predict(lm2), col = 2)
  OUTPUT
} # end of function 





##############################################################################################
##############################################################################################
# B) Computing the topological overlap matrix 
##############################################################################################
##############################################################################################



# ===================================================
#The function TOMdist computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
#
#  ************* Removed: use 1-TOMsimilarity(adjMat). ***********************
#
#TOMdist=function(adjMat, useActualMax = FALSE) 
#{
  #diag(adjMat)=0;
  #adjMat[is.na(adjMat)]=0;
  #maxh1=max(as.dist(adjMat) ); minh1=min(as.dist(adjMat) ); 
  #if (maxh1>1 | minh1 < 0 ) 
    #stop(paste("The adjacency matrix contains entries that are larger than 1 or",
               #"smaller than 0: max =",maxh1,", min =",minh1)) 
  #if ( max(c(as.dist(abs(adjMat-t(adjMat)))))>10^(-12) ) 
    #stop("Non-symmetric adjacency matrix. ") 
  #adjMat= (adjMat+ t(adjMat) )/2
  #connectivity=apply(adjMat,2,sum)
  #maxADJconst=1
  #if (useActualMax==TRUE) maxADJconst=max(c(as.dist(adjMat ))) 
  #Dhelp1=matrix(connectivity,ncol=length(connectivity),nrow=length(connectivity))
  #denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjMat); 
  #gc();gc();
  #numTOM=as.dist(adjMat %*% adjMat +adjMat);
  ##TOMmatrix=numTOM/denomTOM
  ## this turns the TOM matrix into a dissimilarity 
  #out1=1-as.matrix(numTOM/denomTOM) 
  #diag(out1)=1 
  ## setting the diagonal to 1 is unconventional (it should be 0)
  ## but it leads to nicer looking TOM plots... 
  #out1
#}

##---------------------------------------------------------------------------
## This is a somewhat modified TOMdist - most checks are left out as they are
## often not necessary.
#
#  ******* This function is not necessary anymore. Left out. ***********
#
#TOMdistNoChecks = function(adjMat, useActualMax = FALSE)
#{
  #diag(adjMat)=0;
  #adjMat[is.na(adjMat)]=0;
  #connectivity=apply(adjMat,2,sum)
  #maxADJconst=1
  #if (useActualMax==TRUE) maxADJconst=max(c(as.dist(adjMat )))
  #Dhelp1 = matrix(connectivity,ncol=length(connectivity),nrow=length(connectivity))
  #denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(maxADJconst-adjMat);
  #rm(Dhelp1);
  #numTOM=as.dist(adjMat %*% adjMat +adjMat);
  ##TOMmatrix=numTOM/denomTOM
  ## this turns the TOM matrix into a dissimilarity 
  #out1=1-as.matrix(numTOM/denomTOM)
  #rm(numTOM); rm(denomTOM);
  #collectGarbage();
  #diag(out1)=1
  ## setting the diagonal to 1 is unconventional (it should be 0)
  ## but it leads to nicer looking TOM plots... 
  #out1
#}

#---------------------------------------------------------------------------
# exact equivalent of TOMdistNoChecks above, but returns similarity.
# This function works with a generalized adjacency that can be signed.
# If the adjacency is signed, returned TOM will be signed as well (use abs(TOM) to get the usual unsigned
# topological overlap)
# If checkDiag and na.rm are turned both off, the function saves a bit of memory overhead.

# ************* this function is replaced by TOMsimilarity that calls compiled code.

#TOMsimilarity = function(adjMat, useActualMax = FALSE, checkDiag = TRUE, na.rm = TRUE)
#{
  #if (checkDiag) diag(adjMat) = 1;
  #if (na.rm) adjMat[is.na(adjMat)]=0;
  #absAdj = abs(adjMat);
  #connectivity=apply(absAdj,2,sum)-1;
  #maxADJconst=1
  #if (useActualMax==TRUE) maxADJconst=max(c(as.dist(absAdj )))
  #Dhelp1 = matrix(connectivity,ncol=length(connectivity),nrow=length(connectivity))
  #denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(maxADJconst-absAdj);
  #rm(Dhelp1);
  #numTOM=as.dist(adjMat %*% adjMat - adjMat);
  ##TOMmatrix=numTOM/denomTOM
  ## this turns the TOM matrix into a dissimilarity 
  #out1=as.matrix(numTOM/denomTOM)
  #rm(numTOM); rm(denomTOM);
  #collectGarbage();
  #diag(out1)=1
  #out1
#}


# ===================================================
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains binary entries...
# This function is explained in Yip and Horvath (2005)
# http://www.genetics.ucla.edu/labs/horvath/GTOM/
GTOMdist = function(adjMat, degree = 1)
{
  maxh1=max(as.dist(adjMat) ); minh1=min(as.dist(adjMat) );
  if (degree!=round(abs(degree))) 
    stop("'degree' must be a positive integer.");
  if (maxh1>1 | minh1 < 0 )
    stop(paste("Entries of the adjacency matrix are not between 0 and 1: max =",
                 maxh1,", min =",minh1))

  if (  max(c(as.dist(abs(adjMat-t(adjMat)))))>0   ) 
    stop("Given adjacency matrix is not symmetric.")

  B <- adjMat;
  if (degree>=2) for (i in 2:degree) 
  {
          diag(B) <- diag(B) + 1;
          B = B %*% adjMat;# this gives the number of paths with length at most degree connecting a pair
  }   
  B <- (B>0);   # this gives the degree-step reachability from a node to another
  diag(B) <- 0;   # exclude each node being its own neighbor
  B <- B %*% B   # this gives the number of common degree-step-neighbor that a pair of nodes share

  Nk <- diag(B);
  B <- B +adjMat;   # numerator
  diag(B) <- 1;
  denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjMat;
  diag(denomTOM) <- 1;
  1 - B/denomTOM   # this turns the TOM matrix into a dissimilarity
}

#=============================================================================================
#
# vectorTOM: calculate TOM of a vector (or a 'small' matrix) with expression
# data. If the number of columns in vect is small (or 1), number of columns in
# datExpr can be large.
#
#============================================================================================

vectorTOM = function(datExpr, vect, subtract1 = FALSE, blockSize = 2000, 
                     corFnc = "cor", corOptions = "use = 'p'", networkType = "unsigned", power = 6,
                     verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);

  intType = charmatch(networkType, .networkTypes)
  if (is.na(intType))
    stop(paste("Unrecognized 'networkType'. Recognized values are", paste(.networkTypes, collapse = ", ")));

  if (is.null(dim(vect)))
  {
     vect = as.matrix(vect) 
     vectIsVector = TRUE;
  } else vectIsVector = FALSE;

  if (nrow(vect)!=nrow(datExpr)) 
    stop("Input error: numbers of samples in 'vect' and 'datExpr' must be the same.");

  if (ncol(vect)>blockSize) 
    stop(paste("Input error: number of columns in 'vect' is too large. ",
               "If you are certain you want to try anyway, increase 'blockSize' to at least",
               "the number of columns in 'vect'."));

  corEval = parse(text = paste(corFnc, "(datExpr, vect ", prepComma(corOptions), ")"));
  corVE = eval(corEval);
  if (intType==1)
  { corVE = abs(corVE);
  } else if (intType==2)
  { corVE = (1+corVE)/2;
  } else if (intType==3)
  { corVE[corVE < 0] = 0;
  } else 
    stop("Unrecognized networkType argument. Recognized values are 'unsigned', 'signed', and 'signed hybrid'.");

  corVE = corVE^power;

  subtract1 = as.numeric(subtract1);

  nVect = ncol(vect); nGenes = ncol(datExpr);
  TOM = matrix(NA, nrow = nGenes, ncol = nVect);

  if (verbose > 0) {
     if (verbose > 1) cat(paste(spaces, "Calculating TOM of a set of vectors with genes"));
     pind = initProgInd();
  }
  start = 1; 
  denomArr = array(0, dim = c(2, blockSize, nVect));
  while (start <= nGenes)
  {
    end = min(start + blockSize-1, nGenes); 
    blockInd = c(start:end);
    corEval = parse(text = paste(corFnc, "(datExpr[, blockInd], datExpr ", prepComma(corOptions), ")"));
    corEE = eval(corEval);
    if (intType==1)
    { corEE = abs(corEE);
    } else if (intType==2)
    { corEE = (1+corEE)/2;
    } else if (intType==3)
    { corEE[corEE < 0] = 0;
    } 
    corEE = corEE^power;
    num = corEE %*% corVE -subtract1 * corVE[blockInd, ]
    kV = apply(corVE, 2, sum, na.rm = TRUE) - subtract1
    kE = apply(corEE, 1, sum, na.rm = TRUE) - 1;
    denomArr[1, 1:(end-start+1), ] = matrix(kV, nrow = end-start+1, ncol = nVect, byrow = TRUE);
    denomArr[2, 1:(end-start+1), ] = matrix(kE, nrow = end-start+1, ncol = nVect);
    denom = apply(denomArr[, 1:(end-start+1), ], c(2,3), min) + 1 - corVE[blockInd, ];
    TOM[blockInd, ] = num/denom;
    if (verbose > 0) pind = updateProgInd(end/nGenes, pind);
    start = end + 1;
    collectGarbage();
  }
  if (verbose>0) printFlush(" ");

  TOM;
}

#=============================================================================================
#
# subsetTOM: calculate TOM of a subset of vectors with respect to a full set of vectors.
#
#============================================================================================


subsetTOM = function(datExpr, subset, 
                    corFnc = "cor", corOptions = "use = 'p'", networkType = "unsigned", power = 6,
                    verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);

  if (!is.null(dim(subset)))
    stop("'subset' must be a dimensionless vector.");

  if (is.null(dim(datExpr)))
    stop("'datExpr' must be a matrix or data frame.");
  if (length(dim(datExpr))!=2)
    stop("'datExpr' must be two-dimensional.");

  nGenes = ncol(datExpr);

  if (is.logical(subset))
    subset = c(1:nGenes)[subset];

  nBlock = length(subset);

  if (any(!is.finite(subset))) stop("Entries of 'subset' must all be finite.");

  if (min(subset) < 1 | max(subset) > nGenes)
    stop(paste("Some entries of 'subset' are out of range.", 
         "\nNote: 'subset' must contain indices of the subset for which the TOM is calculated."));

  intType = charmatch(networkType, .networkTypes)
  if (is.na(intType))
    stop(paste("Unrecognized 'networkType'. Recognized values are", paste(.networkTypes, collapse = ", ")));

  adj = adjacency(datExpr, subset, power = power, type = networkType, corFnc = corFnc, 
                  corOptions = corOptions);

  adj[is.na(adj)] = 0;
  num = t(adj) %*% adj - adj[subset, ];

  k = apply(adj, 2, sum);

  kMat = matrix(k, nBlock, nBlock);

  denom = pmin(kMat, t(kMat)) - adj[subset, ];

  TOM = num/denom;
  diag(TOM) = 1;

  TOM;
}

#---------------------------------------------------------------------
#
# adjacency
#
#---------------------------------------------------------------------
# Computes the adjacency from the expression data: takes cor, transforms it as appropriate and possibly
# adds a sign if requested. No subselection on datExpr is performed.
# A slighly reworked version that assumes one wants the adjacency matrix of data with itself or a
# subset. The data are given only once, and an additional selection index for columns is given.
# Caution: no checking of selectCols validity is performed.
# The probability method is removed as it's not used.
 
adjacency = function(datExpr, selectCols=NULL, type = "unsigned", power = if (type=="distance") 1 else 6,
                     corFnc = "cor", corOptions = "use = 'p'",
                     distFnc = "dist", distOptions = "method = 'euclidean'")
{
  intType = charmatch(type, .adjacencyTypes)
  if (is.na(intType))
    stop(paste("Unrecognized 'type'. Recognized values are", paste(.adjacencyTypes, collapse = ", ")));

  if (intType < 4)
  {
    if (is.null(selectCols))
    {
      corExpr = parse(text = paste(corFnc, "(datExpr ", prepComma(corOptions), ")"));
      # cor_mat = cor(datExpr, use = "p");
      cor_mat = eval(corExpr);
    } else {
      corExpr = parse(text = paste(corFnc, "(datExpr, datExpr[, selectCols] ", prepComma(corOptions), ")"));
      #cor_mat = cor(datExpr, datExpr[, selectCols], use="p");
      cor_mat = eval(corExpr);
    }
  } else {
    if (!is.null(selectCols)) 
      stop("The argument 'selectCols' cannot be used for distance adjacency.");
    corExpr = parse(text = paste(distFnc, "(t(datExpr) ", prepComma(distOptions), ")"));
    # cor_mat = cor(datExpr, use = "p");
    d = eval(corExpr);
    if (any(d<0)) 
      warning("Function WGCNA::adjacency: Distance function returned (some) negative values.");
    cor_mat = 1-as.matrix( (d/max(d, na.rm = TRUE))^2 );
  }

  if (intType==1)
  { cor_mat = abs(cor_mat); 
  } else if (intType==2)
  { cor_mat = (1+cor_mat)/2; 
  } else if (intType==3)
  { cor_mat[cor_mat < 0] = 0; 
  }
  cor_mat^power;
}

# A presumably faster and less memory-intensive version, only for "unsigned" networks.

unsignedAdjacency = function(datExpr, datExpr2 = NULL, power = 6,
                             corFnc = "cor", corOptions = "use = 'p'")
{
  corExpr = parse(text = paste(corFnc, "(datExpr, datExpr2 ", prepComma(corOptions), ")"));
  # abs(cor(datExpr, datExpr2, use="p"))^power;
  abs(eval(corExpr))^power;
}

#####################################################################################################
#####################################################################################################
# C) Defining gene modules using clustering procedures
#####################################################################################################
#####################################################################################################


cutreeStatic = function(dendro, cutHeight = 0.9, minSize = 50)
{
  normalizeLabels(moduleNumber(dendro, cutHeight, minSize));
}

cutreeStaticColor = function(dendro, cutHeight = 0.9, minSize = 50)
{
  labels2colors(normalizeLabels(moduleNumber(dendro, cutHeight, minSize)));
}

 

plotColorUnderTree = function(
   dendro, 
   colors, 
   rowLabels = NULL, 
   rowWidths = NULL,
   rowText = NULL, 
   textPositions = NULL, 
   addTextGuide = TRUE,
   cex.rowLabels = 1,
   cex.rowText = 0.8,
   ...)
{
  colors = as.matrix(colors);
  dimC = dim(colors)

  if (is.null(rowLabels) & (length(dimnames(colors)[[2]])==dimC[2])) 
     rowLabels = colnames(colors);


  sAF = options("stringsAsFactors")
  options(stringsAsFactors = FALSE);
  on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)

  nColorRows = dimC[2];
  if (length(dendro$order) != dimC[1] ) 
    stop("ERROR: length of colors vector not compatible with number of objects in the hierarchical tree.");
  C = colors[dendro$order, , drop = FALSE]; 
  step = 1/(dimC[1]-1);
  #barplot(height=1, col = "white", border=FALSE, space=0, axes=FALSE, ...)
  barplot(height=1, col = "white", border=FALSE, space=0, axes=FALSE)
  charWidth = strwidth("W")/2;
  if (!is.null(rowText))
  {
     if (is.null(textPositions)) textPositions = c(1:nColorRows);
     if (is.logical(textPositions)) textPositions = c(1:nColorRows)[textPositions];
     nTextRows = length(textPositions);
  } else 
     nTextRows = 0;
  nRows = nColorRows + nTextRows;
  ystep = 1/nRows;
  if (is.null(rowWidths)) 
  { 
    rowWidths = rep(ystep, nColorRows + nTextRows)
  } else
  {
    if (length(rowWidths)!=nRows) 
      stop("plotColorUnderTree: Length of 'rowWidths' must equal the total number of rows.")
    rowWidths = rowWidths/sum(rowWidths);
  }

  hasText = rep(0, nColorRows);
  hasText[textPositions] = 1;
  csPosition = cumsum(c(0, hasText[-nColorRows]));
  
  colorRows = c(1:nColorRows) + csPosition;
  rowType = rep(2, nRows);
  rowType[colorRows] = 1;

  physicalTextRow = c(1:nRows)[rowType==2];

  yBottom = c(0, cumsum(rowWidths[nRows:1])) ;  # Has one extra entry but that shouldn't hurt
  yTop = cumsum(rowWidths[nRows:1]) 

  if (!is.null(rowText))
  {
     charHeight = strheight("W");
     rowText = as.matrix(rowText)
     textPos = list();
     textPosY = list();
     textLevs = list();
     for (tr in 1:nTextRows) 
     {
       width1 = rowWidths[ physicalTextRow[tr] ];
       nCharFit = floor(width1/charHeight/cex.rowText/1.2/par("lheight"))
       if (nCharFit==0) stop("Rows are too narrow to fit text. Consider decreasing cex.rowText.");
       set = textPositions[tr];
       #colLevs = sort(unique(colors[, set]));
       #textLevs[[tr]] = rowText[match(colLevs, colors[, set]), tr];
       textLevs[[tr]] = sort(unique(rowText[, tr]));
       nLevs = length(textLevs[[tr]]);
       textPos[[tr]] = rep(0, nLevs);
       orderedText = rowText[dendro$order, tr]
       for (cl in 1:nLevs)
       {
         ind = orderedText == textLevs[[tr]][cl];
         sind = ind[-1];
         ind1 = ind[-length(ind)];
         starts = c(1:length(ind))[!ind1 & sind];
         ends = c(1:length(ind))[ind1 & !sind];
         if (length(starts)==0) starts = 1;
         if (length(ends)==0) ends = length(ind);
         if (ends[1] < starts[1]) starts = c(1, starts);
         if (ends[length(ends)] < starts[length(starts)]) ends = c(ends, length(ind));
         lengths = ends - starts;
         long = which.max(lengths);
         textPos[[tr]][cl] = starts[long];
       }
       yPos = seq(from = 2, to = 2*nCharFit, by=2) / (2*nCharFit+2);
       textPosY[[tr]] = rep(yPos, ceiling(nLevs/nCharFit)+5)[1:nLevs][rank(textPos[[tr]])];
     }
  } 
        
  jIndex = nRows;

  if (is.null(rowLabels)) rowLabels = c(1:nColorRows);
  C[is.na(C)] = "grey"
  for (j in 1:nColorRows)
  {
    jj = jIndex;
    ind = (1:dimC[1]);
    xl = (ind-1.5) * step; xr = (ind-0.5) * step; 
    yb = rep(yBottom[jj], dimC[1]); yt = rep(yTop[jj], dimC[1]);
    if (is.null(dim(C))) {
       rect(xl, yb, xr, yt, col = as.character(C), border = as.character(C));
    } else {
       rect(xl, yb, xr, yt, col = as.character(C[,j]), border = as.character(C[,j]));
    }
    text(rowLabels[j], pos=2, x= -charWidth -0.5*step, y= (yBottom[jj] + yTop[jj])/2, 
         cex=cex.rowLabels, xpd = TRUE);
    textRow = match(j, textPositions);
    if (is.finite(textRow))
    {
      jIndex = jIndex - 1;
      xt = (textPos[[textRow]] - 1) * step;
      #printFlush(spaste("jIndex: ", jIndex, ", yBottom: ", yBottom[jIndex],
      #                  ", yTop: ", yTop[jIndex], ", min(textPosY): ", min(textPosY[[textRow]]),
      #                  ", max(textPosY): ", max(textPosY[[textRow]])));
      yt = yBottom[jIndex] + (yTop[jIndex]-yBottom[jIndex]) * textPosY[[textRow]];
      nt = length(textLevs[[textRow]]);
      # Add guide lines
      if (addTextGuide)
        for (l in 1:nt) lines(c(xt[l], xt[l]), c(yt[l], yTop[jIndex]), col = "grey70", lty = 3);
      text(textLevs[[textRow]], x = xt, y = yt, adj = c(0, 0.5), xpd = TRUE, cex = cex.rowText)
      # printFlush("ok");
    }
    jIndex = jIndex - 1;
  }
  for (j in 0:(nColorRows + nTextRows)) lines(x=c(0,1), y=c(yBottom[j+1], yBottom[j+1]));
}

#========================================================================================================
# This function can be used to create an average linkage hierarchical
# clustering tree
# or the microarray samples. The rows of datExpr correspond to the samples and
# the columns to the genes
# You can optionally input a quantitative microarray sample trait.

plotClusterTreeSamples=function(datExpr, y = NULL, traitLabels = NULL, yLabels = NULL, 
         main = if (is.null(y)) "Sample dendrogram" else "Sample dendrogram and trait indicator",
         setLayout = TRUE, autoColorHeight = TRUE, colorHeight = 0.3,
         dendroLabels = NULL,
         addGuide = FALSE, guideAll = TRUE, guideCount = NULL,
         guideHang = 0.20, cex.traitLabels = 0.8,
         cex.dendroLabels = 0.9, marAll = c(1,5,3,1),  saveMar = TRUE,
         abHeight = NULL, abCol = "red", ...) 
{
  dendro = flashClust( dist( datExpr  ), method="average" )
  if (is.null(y) ) 
  {
    oldMar = par("mar");
    par(mar = marAll);
    plot(dendro, main=main, sub="", xlab = "", labels = dendroLabels, cex = cex.dendroLabels)
    if (saveMar) par(oldMar);
  } else {
    if (is.null(traitLabels)) traitLabels = names(as.data.frame(y));
    y = as.matrix(y);
    if (!is.numeric(y) ) 
    {
       warning(paste("The microarray sample trait y will be transformed to numeric."));
       dimy = dim(y)
       y=as.numeric(y)
       dim(y) = dimy;
    } # end of if (!is.numeric(y) )
    if (  nrow(as.matrix(datExpr)) != nrow(y) ) 
      stop(paste("Input Error: dim(as.matrix(datExpr))[[1]] != length(y)\n", 
                 "  In plain English: The number of microarray sample arrays does not match the number",
                 "of samples for the trait.\n",
                 "   Hint: Make sure rows of 'datExpr' (and 'y', if it is a matrix) correspond to samples."))

    if (is.integer(y))
    {
      y = y-min(0, min(y, na.rm = TRUE)) + 1;
    } else {
      y = (y>=median(y, na.rm=T)) + 1;
    }
    plotDendroAndColors(dendro, colors = y, groupLabels = traitLabels, rowText = yLabels, 
                        setLayout = setLayout, 
                        autoColorHeight = autoColorHeight, colorHeight = colorHeight,
                        addGuide = addGuide, guideAll = guideAll, guideCount = guideCount, 
                        guideHang = guideHang, cex.colorLabels = cex.traitLabels,
                        cex.dendroLabels = cex.dendroLabels, marAll = marAll, 
                        saveMar = saveMar, abHeight = abHeight, abCol = abCol,
                        main = main,
                        ...);
  }
}# end of function PlotClusterTreeSamples

# ===================================================
# The function TOMplot creates a TOM plot
# Inputs:  distance measure, hierarchical (hclust) object, color label=colors

TOMplot = function(dissim, dendro, Colors=NULL, ColorsLeft = Colors, terrainColors=FALSE, 
                   setLayout = TRUE, ...) 
{
  if ( is.null(Colors) ) Colors=rep("white", dim(as.matrix(dissim))[[1]] )
  if ( is.null(ColorsLeft)) ColorsLeft = Colors;
  nNodes=length(Colors)
  if (nNodes<2) {
     warning("You have only 1 or 2 genes in TOMplot. No plot will be produced")
  } else {
     if (nNodes != length(ColorsLeft)) 
       stop("ERROR: number of (top) color labels does not equal number of left color labels")
     if (nNodes != dim(dissim)[[1]] ) 
       stop(paste("ERROR: number of color labels does not equal number of nodes in dissim.\n",
                  "     nNodes != dim(dissim)[[1]] "))
     labeltree = as.character(Colors)
     labelrow  = as.character(ColorsLeft)
     #labelrow[dendro$order[length(labeltree):1]]=labelrow[dendro$order]
     options(expressions = 10000)
     dendro$height = (dendro$height - min(dendro$height))/(1.15 *
                                     (max(dendro$height)-min(dendro$height)))
     if (terrainColors) {
       .heatmap(as.matrix(dissim), Rowv=as.dendrogram(dendro, hang = 0.1), 
                Colv= as.dendrogram(dendro, hang = 0.1), 
                scale="none", revC=T, ColSideColors=as.character(labeltree),
                RowSideColors=as.character(labelrow), labRow=F, labCol=F, 
                col = terrain.colors(100), setLayout = setLayout, ...) 
     } else {
       .heatmap(as.matrix(dissim), Rowv=as.dendrogram(dendro, hang = 0.1), 
                Colv= as.dendrogram(dendro, hang = 0.1), 
               scale="none",revC=T, ColSideColors=as.character(labeltree),
               RowSideColors=as.character(labelrow), labRow=F, labCol=FALSE, setLayout = setLayout,
               ...)
     } #end of if
  }
} #end of function


plotNetworkHeatmap = function(datExpr,  plotGenes, useTOM = TRUE, power = 6 , 
                              networkType = "unsigned", main = "Heatmap of the network") 
{
  match1=match( plotGenes ,names(data.frame(datExpr)) )
  match1=match1[ !is.na(match1)]
  nGenes=length(match1)
  if (  sum( !is.na(match1) )  != length(plotGenes) ) 
  {
    printFlush(paste("Warning: Not all gene names were recognized.", 
                     "Only the following genes were recognized. "));
    printFlush(paste("   ", names(data.frame(datExpr))[match1], collapse = ", " ))
  }
  if (nGenes< 3 ) 
  { 
    warning(paste("Since you have fewer than 3 genes, the network will not be visualized.\n",
                  "   Hint: please input more genes.")); plot(1,1)
  } else {
    datErest=datExpr[, match1 ]
    ADJ1 = adjacency(datErest, power = power, type = networkType)
    if (useTOM) {
       diss1= 1-TOMsimilarity(ADJ1)   
    } else {
       diss1 = 1-ADJ1;
    }
    diag(diss1)=NA
    hier1=flashClust(as.dist(diss1), method="average" )
    colors1=rep("white", nGenes)
    labeltree = names(data.frame(datErest))
    labelrow  = names(data.frame(datErest))
    labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
    options(expressions = 10000)
    heatmap(as.matrix(diss1),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none", revC=T, 
            labRow= labeltree, labCol= labeltree,main=main)
  } # end of if (nGenes> 2 )
} # end of function

#####################################################################################################
#####################################################################################################
# E) Relating a measure of gene significance to the modules 
#####################################################################################################
#####################################################################################################

# ===================================================
# The function ModuleEnrichment1 creates a bar plot that shows whether modules are enriched with
# significant genes.
# More specifically, it reports the mean gene significance for each module.
# The gene significance can be a binary variable or a quantitative variable.
# It also plots the 95% confidence interval of the mean (CI=mean +/- 1.96* standard error).
# It also reports a Kruskal Wallis P-value.

plotModuleSignificance = function(geneSignificance, colors, boxplot = FALSE, 
                                  main = "Gene significance across modules,",
                                  ylab = "Gene Significance", ...)
{
  if (length(geneSignificance) != length(colors) ) 
    stop("Error: 'geneSignificance' and 'colors' do not have the same lengths")
  no.colors=length(names(table(colors) ))
  if (no.colors==1) pp=NA
  if (no.colors>1) 
  {
    pp=try(kruskal.test(geneSignificance,factor(colors))$p.value)
    if (class(pp)=='try-error') pp=NA
  }
  title = paste(main," p-value=", signif(pp,2), sep = "")
  if (boxplot != TRUE) {
    means1=as.vector(tapply(geneSignificance,colors,mean, na.rm = TRUE));
    se1= as.vector(tapply(geneSignificance,colors,stdErr))
    # par(mfrow=c(1,1))
    barplot(means1, names.arg=names(table(colors) ),col= names(table(colors) ) ,ylab=ylab, 
            main = title, ...)
    addErrorBars(as.vector(means1), as.vector(1.96*se1), two.side=TRUE)
  } else {
    boxplot(split(geneSignificance,colors),notch=T,varwidth=T, col= names(table(colors) ),ylab=ylab,
            main = title, ...)
  }
} # end of function

#####################################################################################################
#####################################################################################################
# F) Carrying out a within module analysis (computing intramodular connectivity etc) 
#####################################################################################################
#####################################################################################################

# ===================================================
#The function DegreeInOut computes for each gene 
#a) the total number of connections, 
#b) the number of connections with genes within its module, 
#c) the number of connections with genes outside its module
# When scaleByMax=TRUE, the within module connectivities are scaled to 1, i.e. the max(K.Within)=1 for each module

intramodularConnectivity = function(adjMat, colors, scaleByMax = FALSE) 
{
  if (nrow(adjMat)!=ncol(adjMat)) 
    stop("'adjMat' is not a square matrix.");
  if (nrow(adjMat)!=length(colors)) 
    stop("Dimensions of 'adjMat' and length of 'colors' differ.");
  nNodes=length(colors)
  colorslevels=levels(factor(colors))
  nLevels=length(colorslevels)
  kWithin=rep(-666,nNodes )
  diag(adjMat)=0
  for (i in c(1:nLevels) ) 
  {
    rest1=colors==colorslevels[i];
    if (sum(rest1) <3 ) { 
       kWithin[rest1]=0 
    } else {
       kWithin[rest1]=apply(adjMat[rest1,rest1], 2, sum, na.rm = TRUE)
       if (scaleByMax) kWithin[rest1]=kWithin[rest1]/max(kWithin[rest1])
    }
  }
  kTotal= apply(adjMat, 2, sum, na.rm = TRUE) 
  kOut=kTotal-kWithin
  if (scaleByMax) kOut=NA
  kDiff=kWithin-kOut
  data.frame(kTotal,kWithin,kOut,kDiff)
}


nPresent = function(x) 
{
  sum(!is.na(x))
}

checkAdjMat = function(adjMat, min = 0, max = 1)
{
  dim = dim(adjMat)
  if (is.null(dim) || length(dim)!=2 )
    stop("adjacency is not two-dimensional");
  if (!is.numeric(adjMat))
    stop("adjacency is not numeric");
  if (dim[1]!=dim[2])
    stop("adjacency is not square");
  if (max(abs(adjMat - t(adjMat)), na.rm = TRUE) > 1e-12)
    stop("adjacency is not symmetric");
  if (min(adjMat, na.rm = TRUE) < min || max(adjMat, na.rm = TRUE) > max)
    stop("some entries are not between", min, "and", max)
}
  


#####################################################################################################
#####################################################################################################
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
#####################################################################################################
#####################################################################################################


# The function signedKME computes the module eigengene based connectivity.
# Input: datExpr= a possibly very large gene expression data set where the rows
# correspond to samples and the columns represent genes
# datME=data frame of module eigengenes (columns correspond to module eigengenes or MEs)
# A module eigengene based connectivity KME value will be computed if the gene has 
# a non missing expression value in at least minNSamples arrays.
# Output a data frame where columns are the KME values 
# corresponding to different modules.
# By splitting the expression data into different blocks, the function can deal with 
# tens of thousands of gene expression data. 
# If there are many eigengenes (say hundreds) consider decreasing the block size.

signedKME = function(datExpr, datME, outputColumnName="kME",
                     corFnc = "cor", corOptions = "use = 'p'") 
{
  datExpr=data.frame(datExpr)
  datME=data.frame(datME)
  output=list()
  if (dim(as.matrix(datME))[[1]] != dim(as.matrix(datExpr))[[1]] ) 
     stop("Number of samples (rows) in 'datExpr' and 'datME' must be the same.")
  varianceZeroIndicatordatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm = TRUE))==0
  varianceZeroIndicatordatME =as.vector(apply(as.matrix(datME),2,var, na.rm = TRUE))==0
  if (sum(varianceZeroIndicatordatExpr,na.rm=T)>0 ) 
    warning("Some genes are constant. Hint: consider removing constant columns from datExpr." )
  if (sum(varianceZeroIndicatordatME,na.rm=T)>0 ) 
    warning(paste("Some module eigengenes are constant, which is suspicious.\n",
            "    Hint: consider removing constant columns from datME." ))
  no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
  if (min(no.presentdatExpr)<..minNSamples ) 
    warning(paste("Some gene expressions have fewer than 4 observations.\n",
            "    Hint: consider removing genes with too many missing values or collect more arrays."))

  #output=data.frame(cor(datExpr, datME, use="p"))
  corExpr = parse(text = paste("data.frame(", corFnc, "(datExpr, datME ", prepComma(corOptions), "))" ));
  output = eval(corExpr);

  output[no.presentdatExpr<..minNSamples, ]=NA
  names(output)=paste(outputColumnName, substring(names(datME), first=3, last=100), sep="")  
  dimnames(output)[[1]] = names(datExpr) 
  output
} # end of function signedKME
 
 


# ===================================================
# The function clusterCoef computes the cluster coefficients.
# Input is an adjacency matrix 

clusterCoef=function(adjMat) 
{
  checkAdjMat(adjMat);
  diag(adjMat)=0
  nNodes=dim(adjMat)[[1]]
  computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
  nolinksNeighbors <- c(rep(-666,nNodes))
  total.edge <- c(rep(-666,nNodes))
  maxh1=max(as.dist(adjMat) ); minh1=min(as.dist(adjMat) ); 
  if (maxh1>1 | minh1 < 0 ) 
    stop(paste("The adjacency matrix contains entries that are larger than 1 or smaller than 0: max =",
                maxh1,", min =",minh1))
  nolinksNeighbors <- apply(adjMat, 1, computeLinksInNeighbors, imatrix=adjMat)
  plainsum  <- apply(adjMat, 1, sum)
  squaresum <- apply(adjMat^2, 1, sum)
  total.edge = plainsum^2 - squaresum
  CChelp=rep(-666, nNodes)
  CChelp=ifelse(total.edge==0,0, nolinksNeighbors/total.edge)
  CChelp
} # end of function



# ===================================================
# The function addErrorBars  is used to create error bars in a barplot
# usage: addErrorBars(as.vector(means), as.vector(stderrs), two.side=F)
addErrorBars<-function(means, errors, two.side=FALSE)
{
 if(!is.numeric(means)) {
      stop("All arguments must be numeric")}

 if(is.null(dim(means)) || length(dim(means))==1){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(means)-1)))) 
 }else{
    if (length(dim(means))==2){
      xval<-cumsum(array(c(1,rep(0,dim(means)[1]-1)),
dim=c(1,length(means))))+0:(length(means)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") }
 }
 MW<-0.25*(max(xval)/length(xval)) 
 ERR1<-means+errors
 ERR2<-means-errors
 for(i in 1:length(means)){
    segments(xval[i],means[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
    if(two.side){
      segments(xval[i],means[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
 } 
} 

# ===================================================
# this function computes the standard error
stdErr <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }

# ===================================================
# The following two functions are for displaying the pair-wise correlation in a panel when using the command "pairs()"
# Typically, we use "pairs(DATA, upper.panel=panel.smooth, lower.panel=.panel.cor, diag.panel=panel.hist)" to
# put the correlation coefficients on the lower panel.

.panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

# ===================================================
# This function is used in "pairs()" function. The problem of the original  panel.cor is that 
# when the correlation coefficient is very small, the lower panel will have a large font 
# instead of a mini-font in a saved .ps file. This new function uses a format for corr=0.2 
# when corr<0.2, but it still reports the original value of corr, with a minimum format.

.panel.cor=function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    txt1=txt
    r1=r
    if (r<0.2) {
        r1=0.2
        txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
        txt1 <- paste(prefix, txt1, sep="")
        }
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt1)
    cex = cex * r1
    r <- round(r, digits)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    text(0.5, 0.5, txt, cex=cex)
}

# ===================================================
# This function collects garbage
# collect_garbage=function(){collectGarbage()}

 
#---------------------------------------------------------------------------------------------------------
# This function plots a barplot with all colors given. If Colors are not given, GlobalStandardColors are
# used, i.e. if you want to see the GlobalStandardColors, just call this function without parameters.
displayColors = function(colors = NULL)
{
  if (is.null(colors)) colors = standardColors();
  barplot(rep(1, length(colors)), col = colors, border = colors);
}


###############################################################################
# I) Functions for merging modules based on a high correlation of the module eigengenes
###############################################################################

#---------------------------------------------------------------------------------------------
#
# dynamicMergeCut
#
#---------------------------------------------------------------------------------------------

dynamicMergeCut = function(n, mergeCor=.9, Zquantile=2.35) 
{
  if (mergeCor>1 | mergeCor<0 ) stop("'mergeCor' must be between 0 and 1.")
  if (mergeCor==1) 
  { 
    printFlush("dynamicMergeCut: given mergeCor=1 will be set to .999."); 
    mergeCor=.999
  }
  if (n<4 ) 
  {
    printFlush(paste("Warning in function dynamicMergeCut: too few observations for the dynamic",
                "assignment of the merge threshold.\n    Will set the threshold to .35"));
    mergethreshold=.35
  } else {
    # Fisher transform of the true merge correlation
    FishermergeCor=.5*log((1+mergeCor)/(1-mergeCor))
    E=exp(2*( FishermergeCor -Zquantile/sqrt(n-3)))
    LowerBoundCIcor=(E-1)/(E+1)
    mergethreshold=1- LowerBoundCIcor
  }
  if (mergethreshold>1) 1 else mergethreshold
}# end of function dynamicMergeCut 



#======================================================================================================
#
# print.flush
#
# =====================================================================================================

#print.flush = function(...)
#{
#   printFlush(...);
#}


##############################################################################################
# I) GENERAL STATISTICAL FUNCTIONS
##############################################################################################

verboseScatterplot = function(x, y, 
                             sample = NULL,
                             corFnc = "cor", corOptions = "use = 'p'",
                             main ="", xlab = NA, ylab = NA, cex=1, cex.axis = 1.5,
                             cex.lab = 1.5, cex.main = 1.5, abline = FALSE, 
                             abline.color = 1, abline.lty = 1,
                             corLabel = corFnc, 
                             displayAsZero = 1e-5,
                             col = 1, bg = 0, 
                             lmFnc = lm,
                             ...) 
{
  if ( is.na(xlab) ) xlab= as.character(match.call(expand.dots = FALSE)$x)
  if ( is.na(ylab) ) ylab= as.character(match.call(expand.dots = FALSE)$y)
  x= as.numeric(as.character(x))
  y= as.numeric(as.character(y))
  corExpr = parse(text = paste(corFnc, "(x, y ", prepComma(corOptions), ")"));
  #cor=signif(cor(x,y,use="p",method=correlationmethod),2)
  cor=signif(eval(corExpr),2)
  if (abs(cor) < displayAsZero) cor = 0;
  corp = signif(corPvalueStudent(cor, sum(is.finite(x) & is.finite(y))), 2);
  #corpExpr = parse(text = paste("cor.test(x, y, ", corOptions, ")"));
  #corp=signif(cor.test(x,y,use="p",method=correlationmethod)$p.value,2)
  #corp=signif(eval(corpExpr)$p.value,2)
  if (corp<10^(-200) ) corp="<1e-200" else corp = paste("=", corp, sep="");
  if (!is.na(corLabel))
  {
     mainX = paste(main, " ", corLabel, "=", cor, ", p",corp, sep="");
  } else
     mainX = main;

  if (!is.null(sample))
  {
    if (length(sample) == 1)
    {
      sample = sample(length(x), sample)
    } 
    if (length(col)<length(x)) col = rep(col, ceiling(length(x)/length(col)));
    if (length(bg )<length(x))  bg = rep(bg,  ceiling(length(x)/length(bg)));
    plot(x[sample], y[sample], main=mainX, xlab=xlab, ylab=ylab, cex=cex, 
         cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main, col = col[sample], bg = bg[sample], ...)
  } else {
    plot(x, y, main=mainX, xlab=xlab, ylab=ylab, cex=cex, 
         cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main, col = col, bg = bg, ...)
  }
  if (abline)
  {
    lmFnc = match.fun(lmFnc);
    fit = lmFnc(y~x);
    abline(reg = fit, col = abline.color, lty = abline.lty);
  }
  invisible(sample);
}

verboseBoxplot = function(x, g,
                          main ="", xlab = NA, ylab = NA, cex=1, cex.axis = 1.5,
                          cex.lab = 1.5, cex.main = 1.5, notch = TRUE, varwidth = TRUE, ...) 
{
  if ( is.na(xlab) ) xlab= as.character(match.call(expand.dots = FALSE)$g)
  #print(xlab1)
  if ( is.na(ylab) ) ylab= as.character( match.call(expand.dots = FALSE)$x)
  #print(ylab1)
  p1 = signif(kruskal.test(x, factor(g) )$p.value,2)
  #if (p1< 5.0*10^(-22) ) p1="< 5e-22"
  boxplot(x~factor(g), notch = notch, varwidth = varwidth,
          main=paste(main,"p =",p1 ),
          xlab=xlab, ylab=ylab, cex=cex, cex.axis=cex.axis,cex.lab=cex.lab, cex.main=cex.main, ...)
}


verboseBarplot = function (x, g,  main = "",
    xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5, cex.lab = 1.5,
    cex.main = 1.5, color="grey", numberStandardErrors=1,
    KruskalTest=TRUE,  AnovaTest=FALSE, two.sided=TRUE, 
    horiz = FALSE, ...) 
{
  stderr1 = function(x){ sqrt( var(x,na.rm=TRUE)/sum(!is.na(x))   ) }
  SE=tapply(x, factor(g), stderr1 )
  err.bp=function(dd,error,two.sided=FALSE,numberStandardErrors, horiz = FALSE)
  {
      if(!is.numeric(dd)) {
            stop("All arguments must be numeric")}
      if(is.vector(dd)){
          xval=(cumsum(c(0.7,rep(1.2,length(dd)-1))))
      }else{
          if (is.matrix(dd)){
            xval=cumsum(array(c(1,rep(0,dim(dd)[1]-1)),
                                dim=c(1,length(dd))))+0:(length(dd)-1)+.5
          }else{
              stop("First argument must either be a vector or a matrix") }
      }
      MW=0.25*(max(xval)/length(xval))
      NoStandardErrors=1
      ERR1=dd+numberStandardErrors*error
      ERR2=dd- numberStandardErrors*error
      if (horiz)
      {
        for(i in 1:length(dd))
        {
          segments(dd[i], xval[i], ERR1[i], xval[i])
          segments(ERR1[i], xval[i]-MW, ERR1[i], xval[i]+MW)
          if(two.sided)
          {
            segments(dd[i], xval[i], ERR2[i], xval[i])
            segments(ERR2[i], xval[i]-MW, ERR2[i], xval[i]+MW)
          }
        }
      } else {
        for(i in 1:length(dd))
        {
          segments(xval[i],dd[i],xval[i],ERR1[i])
          segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
          if(two.sided)
          {
            segments(xval[i],dd[i],xval[i],ERR2[i])
            segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
          }
        }
      }
  } # end of function err.bp
  if ( is.na(ylab) ) ylab= as.character(match.call(expand.dots = FALSE)$x)
  if ( is.na(xlab) ) xlab= as.character( match.call(expand.dots = FALSE)$g)

  Means1= tapply(x, factor(g), mean, na.rm=TRUE)
  if (length(unique(x)) > 2)
  {
    p1=signif(kruskal.test(x~factor(g) )$p.value,2)
    if (AnovaTest)  p1=signif(  anova(lm(x~factor(g)))$Pr[[1]]    ,2)
  } else {
    p1 = tryCatch(signif(fisher.test(x, g, alternative = "two.sided")$p.value, 2), error = function(e) {NA})
  }
  if ( AnovaTest | KruskalTest)  main=paste(main, "p =", p1 )
  ret = barplot(Means1, main=main,col=color, xlab=xlab,ylab=ylab, cex=cex, cex.axis=cex.axis,cex.lab=cex.lab,
                cex.main=cex.main, horiz = horiz, ...)
  abline(h=0)
  if (numberStandardErrors >0) 
  {
      err.bp(as.vector(Means1), as.vector(SE), two.sided=two.sided,
             numberStandardErrors= numberStandardErrors, horiz = horiz)
  }
  attr(ret, "height") = as.vector(Means1);
  attr(ret, "stdErr") = as.vector(SE);
  invisible(ret);
}





#=============================================================================================
#
# Correlation p-value for multiple correlation values
#
#=============================================================================================

corPvalueFisher = function(cor, nSamples, twoSided = TRUE)
{
  if (sum(abs(cor)>1, na.rm = TRUE)>0)
    stop("Some entries in 'cor' are out of normal range -1 to 1.");
  if (twoSided)
  {
     z = abs(0.5 * log((1+cor)/(1-cor)) * sqrt(nSamples-3));
     2 * pnorm(-z);
  } else {
     # return a small p-value for positive correlations
     z = -0.5 * log((1+cor)/(1-cor)) * sqrt(nSamples-3); 
     pnorm(-z);
  }
}

# this function compute an asymptotic p-value for a given correlation (r) and sample size (n) 
# Needs a new name before we commit it to the package.

corPvalueStudent = function(cor, nSamples) 
{
  T=sqrt(nSamples-2) * cor/sqrt(1-cor^2)
  2*pt(abs(T),nSamples-2, lower.tail = FALSE)
}


#########################################################################################

propVarExplained = function(datExpr, colors, MEs, corFnc = "cor", corOptions = "use = 'p'")
{
  fc = as.factor(colors);
  mods = levels(fc);
  nMods = nlevels(fc);
  nGenes = ncol(datExpr);
  if (nMods!=ncol(MEs))
    stop(paste("Input error: number of distinct 'colors' differs from\n", 
               "     the number of module eigengenes given in ME."));

  if (ncol(datExpr)!=length(colors))
    stop("Input error: number of probes (columns) in 'datExpr' differs from the length of goven 'colors'.");

  if (nrow(datExpr)!=nrow(MEs))
    stop("Input error: number of observations (rows) in 'datExpr' and 'MEs' differ.");

  PVE = rep(0, nMods);

  col2MEs = match(mods, substring(names(MEs), 3));

  if (sum(is.na(col2MEs))>0)
    stop("Input error: not all given colors could be matched to names of module eigengenes.");

  for (mod in 1:nMods)
  {
    modGenes = c(1:nGenes)[as.character(colors)==mods[mod]];
    corExpr = parse(text = paste(corFnc, "(datExpr[, modGenes], MEs[, col2MEs[mod]]",
                                 prepComma(corOptions), ")"));
    PVE[mod] = mean(as.vector(eval(corExpr)^2));
  }

  names(PVE) = paste("PVE", mods, sep = "");
  PVE
}
 

#===================================================================================
#
# addGrid
#
#===================================================================================
# This function adds a horizontal grid to a plot 

addGrid = function(linesPerTick = NULL, horiz = TRUE, vert = FALSE, col = "grey30", lty = 3)
{
  box = par("usr");
  if (horiz)
  {
    ticks = par("yaxp");
    nTicks = ticks[3];
    if (is.null(linesPerTick))
    {
       if (nTicks < 6) linesPerTick = 5 else linesPerTick = 2;
    }
    spacing = (ticks[2]-ticks[1])/(linesPerTick*nTicks);
    first = ceiling((box[3] - ticks[1])/spacing);
    last = floor((box[4] - ticks[1])/spacing);
    #print(paste("addGrid: first=", first, ", last =", last, "box = ", paste(signif(box,2), collapse = ", "), 
                #"ticks = ", paste(signif(ticks, 2), collapse = ", "), "spacing =", spacing ));
    for (k in first:last)
      lines(x = box[c(1,2)], y = rep(ticks[1] + spacing * k, 2), 
            col = col, lty = lty);
  }
  if (vert)
  {
    ticks = par("xaxp");
    nTicks = ticks[3];
    if (is.null(linesPerTick))
    {
       if (nTicks < 6) linesPerTick = 5 else linesPerTick = 2;
    }
    spacing = (ticks[2]-ticks[1])/(linesPerTick*ticks[3]);
    first = ceiling((box[1] - ticks[1])/spacing);
    last = floor((box[2] - ticks[1])/spacing);
    #print(paste("addGrid: first=", first, ", last =", last, "box = ", paste(signif(box,2), collapse = ", "), 
    #            "ticks = ", paste(signif(ticks, 2), collapse = ", "), "spacing =", spacing ));
    for (l in first:last)
      lines(x = rep(ticks[1] + spacing * l, 2), y = box[c(3,4)],
            col = col, lty = lty);
  }

}

#-----------------------------------------------------------------------------------------------
#
# Add vertical "guide" lines to a dendrogram to facilitate identification of clusters with color bars
#
#-----------------------------------------------------------------------------------------------

addGuideLines = function(dendro, all = FALSE, count = 50, positions = NULL, col = "grey30", lty = 3,
                         hang = 0)
{
  if (all)
  {
    positions = 1:(length(dendro$height)+1);
  } else {
    if (is.null(positions))
    {
      lineSpacing = (length(dendro$height)+1)/count;
      positions = (1:count)* lineSpacing;
    }
  }
  objHeights = rep(0, length(dendro$height+1));
  objHeights[-dendro$merge[dendro$merge[,1]<0,1]] = dendro$height[dendro$merge[,1]<0];
  objHeights[-dendro$merge[dendro$merge[,2]<0,2]] = dendro$height[dendro$merge[,2]<0];
  box = par("usr"); ymin = box[3]; ymax = box[4];
  objHeights = objHeights - hang*(ymax - ymin);
  objHeights[objHeights<ymin] = ymin;
  posHeights = pmin(objHeights[dendro$order][floor(positions)], objHeights[dendro$order][ceiling(positions)]);
  for (line in 1:length(positions)) # The last guide line is superfluous
    lines(x = rep(positions[line], 2), y = c(ymin, posHeights[line]), lty = 3, col = col);
}

#-------------------------------------------------------------------------------------------
#
# nearestNeighborConnectivity
#
#-------------------------------------------------------------------------------------------
# This function takes expression data (rows=samples, colummns=genes)
# and the power exponent used in weighting the
# correlations to get the network adjacency matrix, and returns an array of dimensions
# nGenes * nSets containing the connectivities of each gene in each subset.

nearestNeighborConnectivity = function(datExpr, nNeighbors = 50, power = 6,
                             type = "unsigned", corFnc = "cor", corOptions = "use = 'p'",
                             blockSize = 1000,  
                             sampleLinks = NULL, nLinks = 5000,
                             setSeed = 38457,
                             verbose=1, indent=0)
{
  spaces = indentSpaces(indent);
  nGenes = dim(datExpr)[2];
  nSamples = dim(datExpr)[1];

  if (is.null(sampleLinks))
  {
    sampleLinks = (nGenes > nLinks);
  }

  if (sampleLinks) nLinks = min(nLinks, nGenes) else nLinks = nGenes;
  
  #printFlush(paste("blockSize =", blockSize));
  #printFlush(paste("nGenes =", nGenes));
  #printFlush(paste(".largestBlockSize =", .largestBlockSize));

  if (blockSize * nLinks>.largestBlockSize) blockSize = as.integer(.largestBlockSize/nLinks);

  intNetworkType = charmatch(type, .networkTypes);
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument. Recognized values are (unique abbreviations of)",
               paste(.networkTypes, collapse = ", ")));

  subtract = rep(1, nGenes);
  if (sampleLinks)
  {
    if (verbose > 0) 
      printFlush(paste(spaces, "nearestNeighborConnectivity: selecting sample pool of size",
                       nLinks, ".."))
    sd = apply(datExpr, 2, sd, na.rm = TRUE);
    order = order(-sd);
    saved = FALSE;
    if (exists(".Random.seed")) 
    {
      saved = TRUE;
      savedSeed = .Random.seed
      if (is.numeric(setSeed)) set.seed(setSeed);
    }
    samplePool = order[sample(x = nGenes, size = nLinks)]
    if (saved) .Random.seed <<- savedSeed;
    poolExpr = datExpr[, samplePool];
    subtract[-samplePool] = 0;
  } 
      
  if (verbose>0) 
  {
     printFlush(paste(spaces, "nearestNeighborConnectivity: received",
                      "dataset with nGenes =", as.character(nGenes)));
     cat(paste(spaces, "..using nNeighbors =", nNeighbors, "and blockSize =", blockSize, "  "));
     pind = initProgInd(trailStr = " done");
  }

  nearestNeighborConn = rep(0, nGenes);

  nBlocks = as.integer((nGenes-1)/blockSize);
  SetRestrConn = NULL;
  start = 1;
  if (sampleLinks)
  {
    corEval = parse(text = paste(corFnc, "(poolExpr, datExpr[, blockIndex] ", prepComma(corOptions), ")"))
  } else {
    corEval = parse(text = paste(corFnc, "(datExpr, datExpr[, blockIndex] ", prepComma(corOptions), ")"))
  }

  while (start <= nGenes)
  {
    end = start + blockSize-1;
    if (end>nGenes) end = nGenes;
    blockIndex = c(start:end);
    #if (verbose>1) printFlush(paste(spaces, "..working on genes", start, "through", end, "of", nGenes))
    c = eval(corEval);
    if (intNetworkType==1)
    { c = abs(c);
    } else if (intNetworkType==2)
    { c = (1+c)/2;
    } else if (intNetworkType==3)
    { c[c < 0] = 0;
    } else stop("Internal error: intNetworkType has wrong value:", intNetworkType, ". Sorry!");
    adj_mat = as.matrix(c^power);
    adj_mat[is.na(adj_mat)] = 0;
    sortedAdj = as.matrix(apply(adj_mat, 2, sort, decreasing = TRUE)[1:(nNeighbors+1), ]);
    nearestNeighborConn[blockIndex] = apply(sortedAdj, 2, sum)-subtract[blockIndex];
    start = end+1;
    if (verbose>0) pind = updateProgInd(end/nGenes, pind);
    collectGarbage();
  }
  if (verbose>0) printFlush(" ");
  nearestNeighborConn;
}


#Try to merge this with the single-set function.
#-------------------------------------------------------------------------------------------
#
# nearestNeighborConnectivityMS
#
#-------------------------------------------------------------------------------------------
# This function takes expression data (rows=samples, colummns=genes) in the multi-set format
# and the power exponent used in weighting the
# correlations to get the network adjacency matrix, and returns an array of dimensions
# nGenes * nSets containing the connectivities of each gene in each subset.

nearestNeighborConnectivityMS = function(multiExpr, nNeighbors = 50, power=6, 
                               type = "unsigned", corFnc = "cor", corOptions = "use = 'p'",
                               blockSize = 1000,
                               sampleLinks = NULL, nLinks = 5000, setSeed = 36492,
                               verbose=1, indent=0)
{
  spaces = indentSpaces(indent);
  setsize = checkSets(multiExpr);
  nGenes = setsize$nGenes;
  nSamples = setsize$nSamples;
  nSets = setsize$nSets;

  if (is.null(sampleLinks))
  {
    sampleLinks = (nGenes > nLinks);
  }

  if (sampleLinks) nLinks = min(nLinks, nGenes) else nLinks = nGenes;

  #printFlush(paste("blockSize =", blockSize));
  #printFlush(paste("nGenes =", nGenes));
  #printFlush(paste(".largestBlockSize =", .largestBlockSize));

  if (blockSize * nLinks>.largestBlockSize) blockSize = as.integer(.largestBlockSize/nLinks);

  if (length(power)==1)
  {
    power = rep(power, nSets);
  } else if (length(power)!=nSets) 
      stop("Invalid arguments: length of 'power' must equal number sets in 'multiExpr'");

  intNetworkType = charmatch(type, .networkTypes);
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument. Recognized values are (unique abbreviations of)",
               paste(.networkTypes, collapse = ", ")));

  subtract = rep(1, nGenes);
  if (sampleLinks)
  { 
    if (verbose > 0) 
      printFlush(paste(spaces, "nearestNeighborConnectivityMS: selecting sample pool of size",
                       nLinks, ".."))
    sd = apply(multiExpr[[1]]$data, 2, sd, na.rm = TRUE);
    order = order(-sd);
    saved = FALSE;
    if (exists(".Random.seed")) 
    {
      saved = TRUE;
      savedSeed = .Random.seed
      if (is.numeric(setSeed)) set.seed(setSeed);
    }
    samplePool = order[sample(x = nGenes, size = nLinks)]
    if (saved) .Random.seed <<- savedSeed;
    subtract[-samplePool] = 0;
  }

  if (verbose>0) printFlush(paste(spaces, "nearestNeighborConnectivityMS: received", nSets, 
                      "datasets with nGenes =", as.character(nGenes)));
  if (verbose>0) printFlush(paste(spaces, "  Using nNeighbors =", nNeighbors));

  nearestNeighborConn = matrix(nrow = nGenes, ncol = nSets);

  if (sampleLinks)
  {
    corEval = parse(text = paste(corFnc, 
          "(multiExpr[[set]]$data[, samplePool], multiExpr[[set]]$data[, blockIndex] ",
                    prepComma(corOptions), ")"))
  } else {
    corEval = parse(text = paste(corFnc, "(multiExpr[[set]]$data, multiExpr[[set]]$data[, blockIndex] ", 
                                 prepComma(corOptions), ")"))
  }


  for (set in 1:nSets) 
  {
    if (verbose>0) {
       cat(paste(spaces, "  Working on set", set));
       pind = initProgInd(trailStr = " done");
    }
    nBlocks = as.integer((nGenes-1)/blockSize);
    SetRestrConn = NULL;
    start = 1;
    while (start <= nGenes)
    {
      end = start + blockSize-1;
      if (end>nGenes) end = nGenes;
      blockIndex = c(start:end);
      #if (verbose>1) printFlush(paste(spaces, " .. working on genes", start, "through", end, "of", nGenes))
      c = eval(corEval);
      if (intNetworkType==1)
      { c = abs(c);
      } else if (intNetworkType==2)
      { c = (1+c)/2;
      } else if (intNetworkType==3)
      { c[c < 0] = 0;
      } else stop("Internal error: intNetworkType has wrong value:", intNetworkType, ". Sorry!");
      adj_mat = as.matrix(c^power[set]);
      adj_mat[is.na(adj_mat)] = 0;
      sortedAdj = as.matrix(apply(adj_mat, 2, sort, decreasing = TRUE)[1:(nNeighbors+1), ]);
      nearestNeighborConn[blockIndex, set] = apply(sortedAdj, 2, sum)-subtract[blockIndex];
      collectGarbage();
      start = end + 1;
      if (verbose > 0) pind = updateProgInd(end/nGenes, pind);
      collectGarbage();
    }
    if (verbose>0) printFlush(" ");
  }
  nearestNeighborConn;
}

#======================================================================================================
#
# Nifty display of progress.
#
# =====================================================================================================

initProgInd = function( leadStr = "..", trailStr = "", quiet = !interactive())
{
  oldStr = " "; 
  cat(oldStr);
  progInd = list(oldStr = oldStr, leadStr = leadStr, trailStr = trailStr);
  class(progInd) = "progressIndicator";
  updateProgInd(0, progInd, quiet);
}

updateProgInd = function(newFrac, progInd, quiet = !interactive())
{
  if (class(progInd)!="progressIndicator") 
    stop("Parameter progInd is not of class 'progressIndicator'. Use initProgInd() to initialize",
         "it prior to use.");

  if (quiet) 
  {
    progInd$oldStr = paste(progInd$leadStr, as.integer(newFrac*100), "% ", progInd$trailStr, sep = "")
  } else {
    cat(paste(rep("\b", nchar(progInd$oldStr)), collapse=""));
    progInd$oldStr = paste(progInd$leadStr, as.integer(newFrac*100), "% ", progInd$trailStr, sep = "")
    cat(progInd$oldStr);
    if (exists("flush.console")) flush.console();
  }
  progInd;
}

#======================================================================================================
#
# Plot a dendrogram and a set of labels underneath
# 
# =====================================================================================================
#

plotDendroAndColors = function(dendro, colors, groupLabels = NULL, rowText = NULL, 
                               textPositions = NULL,
                               setLayout = TRUE, autoColorHeight = TRUE, colorHeight = 0.2,
                               rowWidths = NULL,
                               dendroLabels = NULL, 
                               addGuide = FALSE, guideAll = FALSE, guideCount = 50, 
                               guideHang = 0.20, addTextGuide = FALSE,
                               cex.colorLabels = 0.8, cex.dendroLabels = 0.9,  
                               cex.rowText = 0.8, marAll = c(1,5,3,1),
                               saveMar = TRUE, 
                               abHeight = NULL, abCol = "red", ...)
{
  oldMar = par("mar");
  if (!is.null(dim(colors)))
  {
    nRows = dim(colors)[2];
  } else nRows = 1;
  if (!is.null(rowText)) nRows = nRows + if (is.null(textPositions)) nRows else length(textPositions);
  if (autoColorHeight) colorHeight = 0.2 + 0.3 * (1-exp(-(nRows-1)/6))
  if (setLayout) layout(matrix(c(1:2), 2, 1), heights = c(1-colorHeight, colorHeight));
  par(mar = c(0, marAll[2], marAll[3], marAll[4]));
  plot(dendro, labels = dendroLabels, cex = cex.dendroLabels, ...);
  if (addGuide) 
    addGuideLines(dendro, count = ifelse(guideAll, length(dendro$height)+1, guideCount), hang = guideHang);
  if (!is.null(abHeight)) abline(h=abHeight, col = abCol);
  par(mar = c(marAll[1], marAll[2], 0, marAll[4]));
  plotColorUnderTree(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels, rowText = rowText,
                     textPositions = textPositions, cex.rowText = cex.rowText, rowWidths = rowWidths,
                     addTextGuide = addTextGuide)
  if (saveMar) par(mar = oldMar);
}

####################################################################################################
#
#  Functions included from NetworkScreeningFunctions
#
####################################################################################################

# this function creates pairwise scatter plots between module eigengenes (above the diagonal)
# Below the diagonal are the absolute values of the Pearson correlation coefficients. 
# The diagonal contains histograms of the module eigengene expressions.

plotMEpairs=function(datME, y=NULL, main="Relationship between module eigengenes", clusterMEs=TRUE, ...)
{
  if ( dim(as.matrix(datME))[[2]]==1 & is.null(y) ) 
  {
    hist( datME, ...)
  } else {
    datMEordered=datME
    if (clusterMEs & dim(as.matrix(datME))[[1]] >1 ) 
    {
      dissimME=(1-t(cor(datME, method="p", use="p")))/2
      hclustdatME=flashClust(as.dist(dissimME), method="average" )
      datMEordered=datME[,hclustdatME$order]
    } # end of if
    if ( !is.null(y) ) 
    {
       if ( length(y)  != dim(as.matrix(datMEordered))[[1]] ) 
         stop(paste("The length of the outcome vector 'y' does not match the number of rows of 'datME'.\n",
             "     The columns of datME should correspond to the module eigengenes.\n", 
             "     The rows correspond to the array samples. Hint: consider transposing datME."));
       datMEordered=data.frame(y, datMEordered)
    } # end of if
    pairs( datMEordered,  upper.panel = panel.smooth,     
           lower.panel = .panel.cor , diag.panel=.panel.hist ,main=main, ...)
  } # end if
} # end of function 


#--------------------------------------------------------------------------------------------------
#
# corPredictionSuccess
#
#--------------------------------------------------------------------------------------------------

# The function corPredictionSuccess can be used to determine which method is best for predicting correlations 
# in a new test set. corTestSet should be a vector of correlations in the test set. 
# The parameter topNumber specifies that the top number most positive and the top most negative 
# predicted correlations 
# TopNumber is a vector of integers.
# corPrediction should be a data frame of predictions for the correlations.
# Output a list with the following components:
# meancorTestSetPositive= mean test set correlation among the topNumber of genes 
#    which are predicted to have positive correlations.
# meancorTestSetNegative= mean test set correlation among the topNumber of genes 
#    which are predicted to have negative correlations.
# meancorTestSetOverall=(meancorTestSetPositive-meancorTestSetNegative)/2

corPredictionSuccess=function( corPrediction, corTestSet, topNumber=100 )
{
  nPredictors=dim(as.matrix(corPrediction))[[2]]
  nGenes=dim(as.matrix(corPrediction))[[1]]
  if (length(as.numeric(corTestSet))!=nGenes ) 
     stop("non-compatible dimensions of 'corPrediction' and 'corTestSet'")
  out1=rep(NA, nPredictors)
  meancorTestSetPositive=matrix(NA, ncol=nPredictors, nrow=length(topNumber)  )
  meancorTestSetNegative=matrix(NA, ncol=nPredictors, nrow=length(topNumber)  )
  for (i in c(1:nPredictors) )
  {
    rankpositive=rank(-as.matrix(corPrediction)[,i], ties.method="first")
    ranknegative=rank(as.matrix(corPrediction)[,i], ties.method="first")
    for (j in c(1:length(topNumber) ) ) 
    {
      meancorTestSetPositive[j,i]=mean(corTestSet[rankpositive<= topNumber[j]],na.rm=T)
      meancorTestSetNegative[j,i]= mean(corTestSet[ranknegative<=topNumber[j]],na.rm=T)
    } # end of j loop over topNumber
  } # end of i loop over predictors
  meancorTestSetOverall=data.frame((meancorTestSetPositive-meancorTestSetNegative)/2)
  dimnames(meancorTestSetOverall)[[2]]=names(data.frame(corPrediction)) 
  meancorTestSetOverall=data.frame(topNumber=topNumber, meancorTestSetOverall)
  meancorTestSetPositive=data.frame(meancorTestSetPositive)
  dimnames(meancorTestSetPositive)[[2]]=names(data.frame(corPrediction)) 
  meancorTestSetPositive=data.frame(topNumber=topNumber, meancorTestSetPositive)
  meancorTestSetNegative=data.frame(meancorTestSetNegative)
  dimnames(meancorTestSetNegative)[[2]]=names(data.frame(corPrediction)) 
  meancorTestSetNegative=data.frame(topNumber=topNumber, meancorTestSetNegative)
  datout=list(meancorTestSetOverall=meancorTestSetOverall, meancorTestSetPositive=meancorTestSetPositive, 
              meancorTestSetNegative =meancorTestSetNegative)
  datout
} # end of function corPredictionSuccess



#--------------------------------------------------------------------------------------------------
#
# relativeCorPredictionSuccess
#
#--------------------------------------------------------------------------------------------------

# The function relativeCorPredictionSuccess can be used to test whether a gene screening method 
# is significantly better than a standard method.
# For each gene screening method (column of corPredictionNew) it provides a Kruskal Wallis 
# test p-value for comparison with the vector corPredictionStandard,
# TopNumber is a vector of integers.
# corTestSet should be a vector of correlations in the test set. 
# corPredictionNew should be a data frame of predictions for the 
# correlations.  corPredictionStandard should be the standard prediction (correlation in the training data).
# The function outputs a p-value for the Kruskal test that
# the new correlation prediction methods outperform the standard correlation prediction method.

relativeCorPredictionSuccess=function(corPredictionNew, corPredictionStandard, corTestSet, topNumber=100 )
{
  nPredictors=dim(as.matrix(corPredictionNew))[[2]]
  nGenes=dim(as.matrix(corPredictionNew))[[1]]
  if (length(as.numeric(corTestSet))!=nGenes ) 
     stop("non-compatible dimensions of 'corPrediction' and 'corTestSet'.")
  if (length(as.numeric(corTestSet))!=length(corPredictionStandard) ) 
     stop("non-compatible dimensions of 'corTestSet' and 'corPredictionStandard'.")
  kruskalp=matrix(NA,nrow=length(topNumber), ncol=nPredictors)
  for (i in c(1:nPredictors) )
  {
    rankhighNew=rank(-as.matrix(corPredictionNew)[,i], ties.method="first")
    ranklowNew=rank(as.matrix(corPredictionNew)[,i],ties.method="first")
    for (j in c(1:length(topNumber)) ){
      highCorNew=as.numeric(corTestSet[rankhighNew <= topNumber[j] ])
      lowCorNew=as.numeric(corTestSet[ranklowNew  <= topNumber[j] ])
      highCorStandard=as.numeric(corTestSet[rank(-as.numeric(corPredictionStandard), 
                                                 ties.method="first") <= topNumber[j]])
      lowCorStandard=as.numeric(corTestSet[rank(as.numeric(corPredictionStandard), 
                                                ties.method="first") <= topNumber[j]])
      signedCorNew=c(highCorNew,-lowCorNew)
      signedCorStandard=c(highCorStandard,-lowCorStandard)
      x1=c(signedCorNew,signedCorStandard)
      Grouping=rep(c(2,1), c(length(signedCorNew), length(signedCorStandard)))
      sign1=sign(cor(Grouping,x1, use="p"))
      if (sign1==0) sign1=1
      kruskalp[j,i]=kruskal.test(x=x1, g=Grouping)$p.value*sign1
      #print(names(data.frame(corPredictionNew))[[i]])
      #print(paste("This correlation is positive if the new method is better than the old method" , 
                   # signif(cor(Grouping,x1, use="p"),3)))
    } # end of j loop
  } # end of i loop
  kruskalp[kruskalp<0]=1
  kruskalp=data.frame(kruskalp)
  dimnames(kruskalp)[[2]]= paste(names(data.frame(corPredictionNew)),".kruskalP", sep="")
  kruskalp=data.frame(topNumber=topNumber, kruskalp)
  kruskalp
} # end of function relativeCorPredictionSuccess

#--------------------------------------------------------------------------------------------------
#
# alignExpr
#
#--------------------------------------------------------------------------------------------------

# If y is supplied, it multiplies columns of datExpr by +/-1 to make all correlations with y positive.
# If y is not supplied, the first column of datExpr is used as the reference direction.

alignExpr=function(datExpr, y = NULL) 
{
  if ( !is.null(y) & dim(as.matrix(datExpr))[[1]] != length(y) ) 
    stop("Incompatible number of samples in 'datExpr' and 'y'.")
  if (is.null(y) ) y=as.numeric(datExpr[,1]) 
  sign1=sign(as.numeric(cor(y, datExpr, use="p" )))
  as.data.frame(scale(t(t(datExpr)*sign1)))
} # end of function alignExpr


# this function can be used to rank the values in x. Ties are broken by the method first.
# This function does not appear to be used anywhere in these functions.
#rank1=function(x){
#    rank(x, ties.method="first")
#}

##############################################################################################
#
# Gene expression simulations (functions by P.L.)
#
##############################################################################################

#----------------------------------------------------------------------------
#
# .causalChildren
#
#----------------------------------------------------------------------------
# Note: The returned vector may contain multiple occurences of the same child.

.causalChildren = function(parents, causeMat)
{
  nNodes = dim(causeMat)[[1]];

  # print(paste("Length of parents: ",length(parents)));
  if (length(parents)==0) return(NULL);

  Child_ind = apply(as.matrix(abs(causeMat[, parents])), 1, sum)>0;
  if (sum(Child_ind)>0)
  {
     children = c(1:nNodes)[Child_ind] 
  } else {
     children = NULL;
  }
  children;
}


#----------------------------------------------------------------------------
#
# simulateEigengeneNetwork
#
#----------------------------------------------------------------------------
#
# Given a set of causal anchors, this function creates a network of vectors that should satisfy the
# causal relations encoded in the causal matrix causeMat, i.e. causeMat[j,i] is the causal effect of
# vector i on vector j. 

# The function starts by initializing all vectors to noise given in the noise specification. (The noise
# can be specified for each vector separately.) Then it runs the standard causal network signal
# propagation and returns the resulting vectors.

simulateEigengeneNetwork = function(causeMat, anchorIndex, anchorVectors, noise = 1, verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent);

  if (verbose>0) printFlush(paste(spaces, "Creating seed vectors..."));
  nNodes = dim(causeMat)[[1]];
  nSamples = dim(anchorVectors)[[1]];

  
  if (length(anchorIndex)!=dim(anchorVectors)[[2]])
    stop(paste("Length of anchorIndex must equal the number of vectors in anchorVectors."));

  if (length(noise)==1) noise = rep(noise, nNodes);
  if (length(noise)!=nNodes)
    stop(paste("Length of noise must equal",
            "the number of nodes as given by the dimension of the causeMat matrix."));

  # Initialize all node vectors to noise with given standard deviation

  NodeVectors = matrix(0, nrow = nSamples, ncol = nNodes);
  for (i in 1:nNodes) NodeVectors[,i] = rnorm(n=nSamples, mean=0, sd=noise[i]);

  Levels = rep(0, times = nNodes);

  # Calculate levels for all nodes: start from anchors and go through each successive level of children

  level = 0;
  parents = anchorIndex;
  Children = .causalChildren(parents = parents, causeMat = causeMat);
  if (verbose>1) printFlush(paste(spaces, "..Determining level structure..."));
  while (!is.null(Children))
  {
    # print(paste("level:", level));
    # print(paste("   parents:", parents));
    # print(paste("   Children:", Children));
    level = level + 1;
    if ((verbose>1) & (level/10 == as.integer(level/10))) 
          printFlush(paste(spaces, "  ..Detected level", level));
    #printFlush(paste("Detected level", level));
    Levels[Children] = level;
    parents = Children;
    Children = .causalChildren(parents = parents, causeMat = causeMat);
  }

  HighestLevel = level;

  # Generate the whole network

  if (verbose>1) printFlush(paste(spaces, "..Calculating network..."));
  NodeVectors[,anchorIndex] = NodeVectors[,anchorIndex] + anchorVectors;
  for (level in (1:HighestLevel))
  {
    if ( (verbose>1) & (level/10 == as.integer(level/10)) ) 
      printFlush(paste(spaces, " .Working on level", level));
    #printFlush(paste("Working on level", level));
    LevelChildren = c(1:nNodes)[Levels==level]
    for (child in LevelChildren) 
    {
      LevelParents = c(1:nNodes)[causeMat[child, ]!=0]
      for (parent in LevelParents)
        NodeVectors[, child] = scale(NodeVectors[, child] + causeMat[child, parent]*NodeVectors[,parent]);
    }
  }

  Nodes = list(eigengenes = NodeVectors, causeMat = causeMat, levels = Levels, anchorIndex = anchorIndex);
  Nodes;
}  

#--------------------------------------------------------------------------------------------
#
# simulateModule
#
#--------------------------------------------------------------------------------------------
# The resulting data is normalized.
# Attributes contain the component trueKME giving simulated correlation with module eigengene for
# both module genes and near-module genes. 
# corPower controls how fast the correlation drops with index i in the module; the curve is roughly
# x^{1/corPower} with x<1 and x~0 near the "center", so the higher the power, the faster the curve rises.

simulateModule = function(ME, nGenes, nNearGenes = 0, minCor = 0.3, maxCor = 1, 
                          corPower = 1,
                          signed = FALSE, propNegativeCor = 0.3, geneMeans = NULL, 
                          verbose = 0, indent = 0)
{
    nSamples = length(ME);

    datExpr = matrix(rnorm((nGenes+nNearGenes)*nSamples), nrow = nSamples, 
                            ncol = nGenes+nNearGenes)

    VarME = var(ME)

    # generate the in-module genes
    CorME = maxCor - (c(1:nGenes)/nGenes)^(1/corPower) * (maxCor-minCor);
    noise = sqrt(VarME * (1-CorME^2)/CorME^2);
    sign = rep(1, nGenes);
    if (!signed) 
    {
      negGenes = as.integer(seq(from = 1/propNegativeCor, by = 1/propNegativeCor, 
                                length.out = nGenes * propNegativeCor))
      negGenes = negGenes[negGenes <=nGenes];
      sign[negGenes] = -1;
    }
    for (gene in 1:nGenes)
    {
      datExpr[, gene] = sign[gene] * (ME + rnorm(nSamples, sd = noise[gene]));
    }

    trueKME = CorME;
    # generate the near-module genes

    if (nNearGenes>0) 
    {
      CorME = c(1:nNearGenes)/nNearGenes * minCor;
      noise = sqrt(VarME * (1-CorME^2)/CorME^2);
      sign = rep(1, nNearGenes);
      if (!signed) 
      {
        negGenes = as.integer(seq(from = 1/propNegativeCor, by = 1/propNegativeCor, 
                                  length.out = nNearGenes * propNegativeCor))
        negGenes = negGenes[negGenes <=nNearGenes];
        sign[negGenes] = -1;
      }
      for (gene in 1:nNearGenes)
        datExpr[, nGenes + gene] = ME + sign[gene] * rnorm(nSamples, sd = noise[gene]);
      trueKME = c(trueKME, CorME);
    }

    datExpr = scale(datExpr);
    if (!is.null(geneMeans))
    {
      if (any(is.na(geneMeans)))
        stop("All entries of 'geneMeans' must be finite.");
      if (length(geneMeans)!=nGenes + nNearGenes)
        stop("The lenght of 'geneMeans' must equal nGenes + nNearGenes.");
      datExpr = datExpr + matrix(geneMeans, nSamples, nGenes + nNearGenes, byrow = TRUE);
    }

    attributes(datExpr)$trueKME = trueKME;

    datExpr;

}

#.SimulateModule=function(ME, size,minimumCor=.3) {
#if (size<3) print("WARNING: module size smaller than 3")
#if(minimumCor==0) minimumCor=0.0001;
#maxnoisevariance=var(ME,na.rm=T)*(1/minimumCor^2-1)
#SDvector=sqrt(c(1:size)/size*maxnoisevariance)
#datSignal=suppressWarnings(matrix(c(ME, ME ,-ME),nrow=size ,ncol=length(ME) ,byrow=T))
#datNoise=SDvector* matrix(rnorm(size*length(ME)),nrow=size ,ncol=length(ME))
#datModule=datSignal+datNoise
#t(datModule)
#} # end of function



#--------------------------------------------------------------------------------------------
#
# simulateSmallLayer
#
#--------------------------------------------------------------------------------------------
# Simulates a bunch of small and weakly expressed modules. 

simulateSmallLayer = function(order, nSamples, 
                              minCor = 0.3, maxCor = 0.5, corPower = 1,
                              averageModuleSize, averageExpr, moduleSpacing,
                              verbose = 4, indent = 0)
{
  spaces = indentSpaces(indent);
  nGenes = length(order)
  datExpr = matrix(0, nrow = nSamples, ncol = nGenes);

  maxCorN0 = averageModuleSize;

  if (verbose>0) printFlush(paste(spaces, "simulateSmallLayer: simulating modules with min corr",
          minCor, ", average expression", averageExpr, ", average module size", averageModuleSize, 
          ", inverse density", moduleSpacing));

  index = 0;
  while (index < nGenes)
  {
     ModSize = as.integer(rexp(1, 1/averageModuleSize));
     if (ModSize<3) ModSize = 3;
     if (index + ModSize>nGenes) ModSize = nGenes - index;
     if (ModSize>2)   # Otherwise don't bother :)
     {
       ModuleExpr = rexp(1, 1/averageExpr);
       if (verbose>4) printFlush(paste(spaces, "  Module of size", ModSize, ", expression", ModuleExpr, 
                                  ", min corr", minCor, 
                                  "inserted at index", index+1));
       ME = rnorm(nSamples, sd = ModuleExpr);
       NInModule = as.integer(ModSize*2/3);
       nNearModule = ModSize - NInModule;
       EffMinCor = minCor * maxCor;
       datExpr[, order[(index+1):(index + ModSize)]] = 
           ModuleExpr * simulateModule(ME, NInModule, nNearModule, EffMinCor, maxCor, corPower);
     }
     index = index + ModSize * moduleSpacing;
  }
  datExpr;
}
         
     
#--------------------------------------------------------------------------------------------
#
# simulateDatExpr
#
#--------------------------------------------------------------------------------------------
#
# Caution: the last Mod.Props entry gives the number of "true grey" genes;
# the corresponding minCor entry must be absent (i.e. length(minCor) = length(modProportions)-1

# SubmoduleLayers: layers of small modules with weaker correlation, ordered in the same order as the
# genes in the big modules. Needs average number of genes in a module (exponential distribution),
# average expression strength (exponential density) and inverse density.

# ScatteredModuleLayers: Layers of small modules whose order is random.

simulateDatExpr=function(eigengenes, nGenes, modProportions,
                          minCor = 0.3, maxCor = 1, 
                          corPower = 1, 
                          signed = FALSE, propNegativeCor = 0.3,
                          geneMeans = NULL,
                          backgroundNoise = 0.1, leaveOut = NULL,
			  nSubmoduleLayers = 0, nScatteredModuleLayers = 0, 
                          averageNGenesInSubmodule = 10, averageExprInSubmodule = 0.2, 
                          submoduleSpacing = 2,
                          verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);

  nMods=length(modProportions)-1;

  nSamples = dim(eigengenes)[[1]];

  if (length(minCor)==1) minCor = rep(minCor, nMods);
  if (length(maxCor)==1) maxCor = rep(maxCor, nMods);

  if (length(minCor)!=nMods)
    stop(paste("Input error: minCor is an array of different lentgh than",
                "the length-1 of modProportions array."));

  if (length(maxCor)!=nMods)
    stop(paste("Input error: maxCor is an array of different lentgh than",
                "the length-1 of modProportions array."));

  if (dim(eigengenes)[[2]]!=nMods)
     stop(paste("Input error: Number of seed vectors must equal the",
                "length of modProportions."));

  if (is.null(geneMeans)) geneMeans = rep(0, nGenes);
  if (length(geneMeans)!=nGenes)
    stop("Length of 'geneMeans' must equal 'nGenes'.");
 
  if (any(is.na(geneMeans)))
    stop("All entries of 'geneMeans' must be finite.");
       
  grey = 0;
  moduleLabels = c(1:nMods);

  if(sum(modProportions)>1) stop("Input error: the sum of Mod.Props must be less than 1");
  #if(sum(modProportions[c(1:(length(modProportions)-1))])>=0.5) 
         #print(paste("SimulateExprData: Input warning: the sum of modProportions for proper modules",
                                       #"should ideally be less than 0.5."));

  no.in.modules = as.integer(nGenes*modProportions);
  no.in.proper.modules = no.in.modules[c(1:(length(modProportions)-1))];
  no.near.modules = as.integer((nGenes - sum(no.in.modules)) * 
                         no.in.proper.modules/sum(no.in.proper.modules));

  simulate.module = rep(TRUE, times = nMods);
  if (!is.null(leaveOut)) simulate.module[leaveOut] = FALSE;

  no.in.modules[nMods+1] = nGenes - sum(no.in.proper.modules[simulate.module]) -
                                          sum(no.near.modules[simulate.module]);

  labelOrder = moduleLabels[rank(-modProportions[-length(modProportions)], ties.method = "first")];
  labelOrder = c(labelOrder, grey);

  if (verbose>0) printFlush(paste(spaces, "simulateDatExpr: simulating", nGenes, "genes in",
                        nMods, "modules."));

  if (verbose>1) 
  {
  #  printFlush(paste(spaces, "    Minimum correlation in a module is", minCor, 
  #                            " and its dropoff is characterized by power", corPower));
    printFlush(paste(spaces, "    Simulated labels:", 
                       paste(labelOrder[1:nMods], collapse = ", "), " and ", grey));
    printFlush(paste(spaces, "    Module sizes:", paste(no.in.modules, collapse = ", ")));
    printFlush(paste(spaces, "    near module sizes:", paste(no.near.modules, collapse = ", ")));
    printFlush(paste(spaces, "    Min correaltion:", paste(minCor, collapse = ", ")));
    if (!is.null(leaveOut)) printFlush(paste(spaces, "    _leaving out_ modules", 
                                              paste(labelOrder[leaveOut], collapse = ", ")));
    
  }

  truemodule=rep(grey, nGenes);
  allLabels=rep(grey, nGenes);	# These have the colors for left-out modules as well.
  
  # This matrix contains the simulated expression values (rows are genes, columns samples)
  # Each simulated cluster has a distinct mean expression across the samples

  datExpr = matrix(rnorm(nGenes*nSamples), nrow = nSamples, ncol = nGenes)
  trueKME = rep(NA, nGenes);
  trueKME.whichMod = rep(0, nGenes);

  gene.index = 0;		# Where to put the current gene into datExpr

  for(mod in c(1:nMods)) 
  {
     nModGenes = no.in.modules[mod];
     nNearGenes = no.near.modules[mod];
     if (simulate.module[mod])
     {
       ME = eigengenes[, mod];
       EffMaxCor = maxCor[mod]; 
       EffMinCor = minCor[mod]; 
       range = (gene.index+1):(gene.index+nModGenes+nNearGenes);
       temp = simulateModule(ME, nModGenes, nNearGenes, minCor[mod], maxCor[mod], 
                         corPower, 
                         signed = signed, propNegativeCor = propNegativeCor,
                         geneMeans = NULL,
                         verbose = verbose-2, indent = indent+2);
       datExpr[, range] = temp;
       truemodule[(gene.index+1):(gene.index+nModGenes)] = labelOrder[mod];
       trueKME[range] = attributes(temp)$trueKME;
       trueKME.whichMod[range] = mod;
     } 
     allLabels[(gene.index+1):(gene.index+nModGenes)] = labelOrder[mod];
     gene.index = gene.index + nModGenes + nNearGenes;
  }

  if (nSubmoduleLayers>0) 
  {
    OrderVector = c(1:nGenes)
    for (layer in 1:nSubmoduleLayers)
    {
      if (verbose>1) printFlush(paste(spaces, "Simulating ordereded extra layer", layer)); 
      datExpr = datExpr + simulateSmallLayer(OrderVector, nSamples, minCor[1], 
                                    maxCor[1],
                                    corPower, averageNGenesInSubmodule, 
                                    averageExprInSubmodule, submoduleSpacing,
                                    verbose-1, indent+1);
    }
  }
  if (nScatteredModuleLayers>0) for (layer in 1:nScatteredModuleLayers)
  {
    if (verbose>1) printFlush(paste(spaces, "Simulating unordereded extra layer", layer)); 
    OrderVector = sample(nGenes)
    datExpr = datExpr + simulateSmallLayer(OrderVector, nSamples, minCor[1],
                                    maxCor[1], corPower, 
                                    averageNGenesInSubmodule, 
                                    averageExprInSubmodule, submoduleSpacing, 
                                    verbose = verbose-1, indent = indent+1);
  }
  collectGarbage();
  if (verbose>1) printFlush(paste(spaces, "  Adding background noise with amplitude", backgroundNoise));
  datExpr = datExpr + rnorm(n = nGenes*nSamples, sd = backgroundNoise);
  means = colMeans(datExpr);

  datExpr = datExpr + matrix(geneMeans - means, nSamples, nGenes, byrow = TRUE);

  colnames(datExpr) = spaste("Gene.", c(1:nGenes));
  rownames(datExpr) = spaste("Sample.", c(1:nSamples));

  list(datExpr = datExpr, setLabels = truemodule, allLabels = allLabels, 
       labelOrder = labelOrder, trueKME = trueKME, trueKME.whichMod = trueKME.whichMod)
} # end of function

#--------------------------------------------------------------------------------------
#
# simulateMultiExpr
#
#--------------------------------------------------------------------------------------
# simulate several sets with some of the modules left out. 
# eigengenes are specified in a standard multi-set data format.
# leaveOut must be a matrix of No.Modules x No.Sets of TRUE/FALSE values;
# minCor must be a single number here; modProportions are a single vector, since the proportions should be the
# same for all sets.
# nSamples is a vector specifying the number of samples in each set; this must be compatible with the
# dimensions of the eigengenes.

simulateMultiExpr = function(eigengenes, nGenes, modProportions,
                          minCor = 0.5, maxCor = 1, 
                          corPower = 1, backgroundNoise = 0.1, leaveOut = NULL,
                          signed = FALSE, propNegativeCor = 0.3,
                          geneMeans = NULL,
			  nSubmoduleLayers = 0, nScatteredModuleLayers = 0, 
                          averageNGenesInSubmodule = 10, averageExprInSubmodule = 0.2, 
                          submoduleSpacing = 2,
                          verbose = 1, indent = 0)
{
  MEsize = checkSets(eigengenes);
  nSets = MEsize$nSets;
  nMods = MEsize$nGenes;
  nSamples = MEsize$nSamples;

  nAllSamples = sum(nSamples);

  if (is.null(geneMeans))
  {
     geneMeans = matrix(0, nGenes, nSets);
  } else {
     geneMeans = as.matrix(geneMeans);
     if (nrow(geneMeans)!=nGenes)
     {
       stop("Number of rows (or entries) in 'geneMeans' must equal 'nGenes'.");
     } else if (ncol(geneMeans)==1)
     {
        geneMeans = matrix(geneMeans, nGenes, nSets);
     } else if (ncol(geneMeans)!=nSets)
        stop("Number of columns in geneMeans must either equal the number of sets or be 1.");
  }

  if (any(is.na(geneMeans)))
    stop("All entries of 'geneMeans' must be finite.");
       
  d2 = length(modProportions)-1;
  if (d2 != nMods) stop(paste("Incompatible numbers of modules in 'eigengenes' and 'modProportions'"));
  if (is.null(leaveOut))
  {
    leaveOut = matrix(FALSE, nMods, nSets);
  } else {
    d3 = dim(leaveOut);
    if ( (d3[1] != nMods) | (d3[2] != nSets) ) 
      stop(paste("Incompatible dimensions of 'leaveOut' and set eigengenes."))
  }

  multiExpr = vector(mode="list", length = nSets);
  setLabels = NULL;
  allLabels = NULL;
  labelOrder = NULL;

  for (set in 1:nSets)
  {
    SetEigengenes = scale(eigengenes[[set]]$data);
    setLeaveOut = leaveOut[, set];
    # Convert setLeaveOut from boolean to a list of indices where it's TRUE
    # SetMinCor = rep(minCor, nMods);
    # SetMaxCor = rep(maxCor, nMods);
    SetLO = c(1:nMods)[setLeaveOut];
    setData = simulateDatExpr(SetEigengenes, nGenes, modProportions,
                          minCor = minCor, maxCor = maxCor, 
                          corPower = corPower, 
                          signed = signed, propNegativeCor = propNegativeCor,
                          backgroundNoise = backgroundNoise, leaveOut = SetLO,
			  nSubmoduleLayers = nSubmoduleLayers,
                          nScatteredModuleLayers  = nScatteredModuleLayers , 
                          averageNGenesInSubmodule = averageNGenesInSubmodule, 
                          averageExprInSubmodule = averageExprInSubmodule, 
                          submoduleSpacing = submoduleSpacing,
                          verbose = verbose-1, indent = indent+1);
    multiExpr[[set]] = list(data = setData$datExpr);
    setLabels = cbind(setLabels, setData$setLabels);
    allLabels = cbind(allLabels, setData$allLabels);
    labelOrder = cbind(labelOrder, setData$labelOrder);
  }
  list(multiExpr = multiExpr, setLabels = setLabels, allLabels = allLabels, 
       labelOrder = labelOrder);
} 

#--------------------------------------------------------------------------------------------------
#
# simulateDatExpr5Modules 
#
#--------------------------------------------------------------------------------------------------

simulateDatExpr5Modules = function(
     nGenes=2000, 
     colorLabels=c("turquoise","blue", "brown", "yellow", "green"),
     simulateProportions=c(0.10,0.08, 0.06, 0.04, 0.02),
     MEturquoise,
     MEblue,
     MEbrown,
     MEyellow,
     MEgreen,
     SDnoise=1,   
     backgroundCor=0.3)
{
   nSamples=length(MEturquoise)
   if( length(MEturquoise) != length(MEblue) | length(MEturquoise) != length(MEbrown) | 
       length(MEturquoise) != length(MEyellow) | length(MEturquoise) != length(MEgreen) ) 
     stop("Numbers of samples in module eigengenes (MEs) are not consistent" );
   if ( sum(simulateProportions)>1 ) 
   { 
     stop("Sum of module proportions is larger than 1. Please ensure sum(simulateProportions)<=1. " ); 
     # simulateProportions=rep(1/10,5)
   } 
   modulesizes=round(nGenes*c(simulateProportions, 1-sum(simulateProportions)))
   truemodule=rep(c( as.character(colorLabels),"grey") , modulesizes )
   ModuleEigengenes = data.frame(MEturquoise,MEblue,MEbrown,MEyellow,MEgreen)
   no.MEs=dim(ModuleEigengenes)[[2]]
   # This matrix contains the simulated expression values 
   #(rows are samples, columns genes)
   # it contains some background noise 
   datExpr=matrix(rnorm(nSamples*nGenes,mean=0,sd=SDnoise),nrow=nSamples,ncol=nGenes)

   if (is.logical(backgroundCor)) backgroundCor = 0.3 * as.numeric(backgroundCor);

   if (as.numeric(backgroundCor) > 0)  
   {
     MEbackground=MEturquoise
     datSignal= (matrix(MEbackground,nrow=length(MEturquoise) ,ncol=nGenes,byrow=F))
     datExpr= datExpr+ as.numeric(backgroundCor)*datSignal
   }# end of if backgroundCor

   for (i in c(1:no.MEs) ) 
   {
     restrict1= truemodule== colorLabels[i]
     datModule = simulateModule(ModuleEigengenes[,i] , nGenes = modulesizes[i], corPower = 2.5)
     datExpr[,restrict1]= datModule
   } # end of for loop
   # this is the output of the function
   list(datExpr =datExpr, truemodule =truemodule, datME = ModuleEigengenes ) 
} # end of simulation function


#--------------------------------------------------------------------------------------------------
#
# automaticNetworkScreening
#
#--------------------------------------------------------------------------------------------------


automaticNetworkScreening = function(
       datExpr, y,   
       power=6, 
       networkType="unsigned", 
       detectCutHeight = 0.995,
       minModuleSize = min(20, ncol(as.matrix(datExpr))/2 ), 
       datME=NULL,  
       getQValues = TRUE, ...) 
{
  y = as.numeric(as.character(y))
  if (length(y) != dim(as.matrix(datExpr))[[1]] ) 
    stop("Number of samples in 'y' and 'datExpr' disagree: length(y) != dim(as.matrix(datExpr))[[1]] ")

  nAvailable=apply(as.matrix(!is.na(datExpr)), 2,sum)
  ExprVariance=apply(as.matrix(datExpr),2,var, na.rm=T ) 
  restrictGenes = (nAvailable>=..minNSamples) & (ExprVariance>0)
  numberUsefulGenes=sum(restrictGenes,na.rm=T) 
  if ( numberUsefulGenes<3 ) 
  {
    stop(paste("IMPORTANT: there are not enough useful genes. \n", 
       "    Your input genes have fewer than 4 observations or they are constant.\n",
       "    WGCNA cannot be used for these data. Hint: collect more arrays or input genes that vary."));
    #warning(paste("IMPORTANT: there are not enough useful genes. \n", 
    #   "    Your input genes have fewer than 4 observations or they are constant.\n",
    #   "    WGCNA cannot be used for these data. Hint: collect more arrays or input genes that vary."));
    #output=list(NetworkScreening=data.frame(NS1=rep(NA, dim(as.matrix(datExpr))[[2]] )), 
    #            datME=rep(NA, dim(as.matrix(datExpr))[[1]] ), EigengeneSignificance=NA , AAcriterion=NA)
    #return(output);
  }

  datExprUsefulGenes=as.matrix(datExpr)[,restrictGenes & !is.na(restrictGenes)]
  if (is.null(datME) )
  {
    mergeCutHeight1 = dynamicMergeCut(n= dim(as.matrix(datExprUsefulGenes))[[1]])
    B = blockwiseModules(datExprUsefulGenes, mergeCutHeight = mergeCutHeight1,  
                         TOMType = "none", power = power, networkType=networkType,
                         detectCutHeight = detectCutHeight, minModuleSize = minModuleSize );
    datME=data.frame(B$MEs)
  }

  if (dim(as.matrix(datME))[[1]] != dim(as.matrix(datExpr))[[1]] ) 
     stop(paste("Numbers of samples in 'datME' and 'datExpr' are incompatible:", 
          "dim(as.matrix(datME))[[1]] != dim(as.matrix(datExpr))[[1]]"))

  MMdata=signedKME(datExpr=datExpr, datME=datME, outputColumnName="MM.")
  MMdataPvalue=as.matrix(corPvalueStudent(as.matrix(MMdata), nSamples= dim(as.matrix(datExpr))[[1]]))
  dimnames( MMdataPvalue)[[2]]=paste("Pvalue",names(MMdata), sep=".")

  NS1=networkScreening(y= y,datME=datME, datExpr=datExpr, getQValues = getQValues)
  # here we compute the eigengene significance measures
  ES=data.frame(cor(y, datME, use="p"))

  ESvector = as.vector(as.matrix(ES));
  EScounts = tapply(abs(ESvector),cut(abs(ESvector),seq(from=0,to=1, by=.1)),length )
  EScounts[is.na(EScounts)] = 0;

  rr=max(abs(ES),na.rm=T)
  AAcriterion=sqrt(length(y)-2) * rr/sqrt(1-rr^2)


  ESy=(1+max(abs(ES), na.rm=T))/2
  ES=data.frame(ES, ESy=ESy)
  
  # to avoid dividing by zero, we set correlation that are 1 equal to .9999
  ES.999=as.numeric(as.vector(ES))
  ES.999[!is.na(ES) &  ES>0.9999]=.9999
  ES.pvalue=corPvalueStudent(cor=abs(ES.999), nSamples=sum(!is.na(y) )) 
  ES.pvalue[length(ES.999)]=0
  EigengeneSignificance.pvalue=data.frame(matrix(ES.pvalue, nrow=1)   )
  names(EigengeneSignificance.pvalue)=names(ES)

  datME=data.frame(datME,y=y)
  names(ES)=paste("ES", substr(names(ES),3,100), sep="")
  
  print(signif(ES,2))

  output=list(networkScreening=data.frame(NS1, MMdata, MMdataPvalue), datME=data.frame(datME), 
              eigengeneSignificance=data.frame(ES) , 
              EScounts = EScounts,
              eigengeneSignificance.pvalue=EigengeneSignificance.pvalue, 
              AAcriterion=AAcriterion)
  
  output
} # end of function automaticNetworkScreening


#--------------------------------------------------------------------------------------------------
#
# automaticNetworkScreeningGS
#
#--------------------------------------------------------------------------------------------------

automaticNetworkScreeningGS = function(
         datExpr, GS,   
         power=6, networkType="unsigned",  detectCutHeight = 0.995,
         minModuleSize = min(20, ncol(as.matrix(datExpr))/2 ), datME=NULL) 
{
  if (!is.numeric(GS) ) 
     stop("Gene significance 'GS' is not numeric.")
  if (  dim(as.matrix(datExpr))[[2]] != length(GS) ) 
     stop("length of gene significance variable GS does not equal the number of columns of datExpr.");

  mergeCutHeight1 = dynamicMergeCut(n= dim(as.matrix(datExpr))[[1]])
  nAvailable=apply(as.matrix(!is.na(datExpr)), 2,sum)
  ExprVariance=apply(as.matrix(datExpr),2,var, na.rm=T ) 
  restrictGenes=nAvailable>=4 & ExprVariance>0
  numberUsefulGenes=sum(restrictGenes,na.rm=T) 
  if ( numberUsefulGenes<3 ) 
  {
    stop(paste("IMPORTANT: there are not enough useful genes. \n", 
       "    Your input genes have fewer than 4 observations or they are constant.\n",
       "    WGCNA cannot be used for these data. Hint: collect more arrays or input genes that vary."));
    #output=list(NetworkScreening=data.frame(NS1=rep(NA, dim(as.matrix(datExpr))[[2]]))  , datME=rep(NA,
    #dim(as.matrix(datExpr))[[1]])    , hubGeneSignificance=NA);
  } # end of if 
  datExprUsefulGenes=as.matrix(datExpr)[,restrictGenes & !is.na(restrictGenes)]

  if (is.null(datME) )
  {
     B = blockwiseModules(datExprUsefulGenes, mergeCutHeight = mergeCutHeight1,  
                        TOMType = "none", power = power, networkType = networkType,
                        detectCutHeight = detectCutHeight, minModuleSize= minModuleSize );
     datME = data.frame(B$MEs)
  } #end of if
  MMdata=signedKME(datExpr=datExpr, datME=datME, outputColumnName="MM.")
  MMdataPvalue=as.matrix(corPvalueStudent(as.matrix(MMdata), nSamples= dim(as.matrix(datExpr))[[1]]))
  dimnames( MMdataPvalue)[[2]]=paste("Pvalue",names(MMdata), sep=".")
  
  NS1= networkScreeningGS(datExpr=datExpr, datME=datME,  GS=GS )
  # here we compute the eigengene significance measures
  HGS1=data.frame(as.matrix(t(hubGeneSignificance(MMdata ^3,GS^3)),nrow=1))
  datME=data.frame(datME)
  names(HGS1)=paste("HGS", substr(names(MMdata),4,100), sep="")
  # now we compute the AA criterion
  print(signif(HGS1,2))
  output = list(networkScreening=data.frame(NS1, MMdata, MMdataPvalue), datME=data.frame(datME), 
                hubGeneSignificance=data.frame(HGS1))
  output
} # end of function automaticNetworkScreeningGS


#--------------------------------------------------------------------------------------------
#
#  hubGeneSignificance
#
#--------------------------------------------------------------------------------------------

# The following function computes the hub gene significance as defined in
# in the paper Horvath and Dong. Input a data frame with possibly signed
# module membership measures ( also known as module eigengene based connectivity
#kME. Further it requires a possibly signed gene significance measure.
# GS=0 means that the gene is not significant, high positive or negative values mean
# that it is significant.
# The input to this function can include the sign of the correlation.
hubGeneSignificance=function(datKME, GS ) 
{
  nMEs=dim(as.matrix(datKME))[[2]]
  nGenes= dim(as.matrix(datKME))[[1]]
  if ( length(GS) !=  nGenes ) 
    stop("Numbers of genes in 'datKME' and 'GS' are not compatible. ")
  Kmax=as.numeric(apply(as.matrix(abs(datKME)),2,max, na.rm=T))
  Kmax[Kmax==0]=1
  datKME=scale(datKME, center=F, scale=Kmax)
  sumKsq=as.numeric(apply(as.matrix(datKME^2) , 2, sum, na.rm=T))
  sumKsq[sumKsq==0]=1
  HGS=as.numeric(apply(I(GS)*datKME, 2, sum,na.rm=T))/ sumKsq
  as.numeric(HGS)
} #end of function hubGeneSignificance


#--------------------------------------------------------------------------------------------
#
#  networkScreeningGS
#
#--------------------------------------------------------------------------------------------

networkScreeningGS = function(datExpr , datME, GS ,
           oddPower = 3, 
           blockSize = 1000,
           minimumSampleSize = ..minNSamples,
           addGS=TRUE)
{
  oddPower=as.integer(oddPower)
  if (as.integer(oddPower/2)==oddPower/2 ) {oddPower=oddPower+1}
  nMEs=dim(as.matrix(datME))[[2]]
  nGenes=dim(as.matrix(datExpr))[[2]]
  GS.Weighted=rep(0,nGenes)

  if ( dim(as.matrix(datExpr))[[1]] != dim(as.matrix(datME))[[1]]) 
    stop(paste("Expression data and the module eigengenes have different\n",
               "      numbers of observations (arrays). Specifically:\n",
               "      dim(as.matrix(datExpr))[[1]] != dim(as.matrix(datME))[[1]] "))

  if ( dim(as.matrix(datExpr))[[2]] != length(GS) ) 
    stop(paste("The number of genes in the expression data does not match\n",
           "      the length of the genes significance variable. Specifically:\n",
           "       dim(as.matrix(datExpr))[[2]] != length(GS)   "));

  nAvailable=apply(as.matrix(!is.na(datExpr)), 2,sum)
  ExprVariance=apply(as.matrix(datExpr),2,var, na.rm=T ) 
  restrictGenes=nAvailable>=4 & ExprVariance>0
  numberUsefulGenes=sum(restrictGenes,na.rm=T) 
  if ( numberUsefulGenes<3 ) 
  {
    stop(paste("IMPORTANT: there are fewer than 3 useful genes. \n", 
       "    Violations: either fewer than 4 observations or they are constant.\n",
       "    WGCNA cannot be used for these data. Hint: collect more arrays or input genes that vary."));
    # datout=data.frame(GS.Weighted=rep(NA, dim(as.matrix(datExpr))[[2]]), GS=GS)
  } # end of if 

  nBlocks=as.integer(nMEs/blockSize)
  if (nBlocks>0) for (i in 1:nBlocks) 
  {
    printFlush(paste("block number = ", i))
    index1=c(1:blockSize)+(i-1)* blockSize
    datMEBatch= datME[,index1]
    datKMEBatch=as.matrix(signedKME(datExpr,datMEBatch, outputColumnName="MM."))
    ESBatch=   hubGeneSignificance(datKMEBatch ^oddPower,GS^oddPower)
    # the following omits the diagonal when datME=datExpr
    if (nGenes==nMEs) {diag(datKMEBatch[index1,])=0
      # missing values will not be used 
      datKMEBatch[is.na(datKMEBatch)]=0
      ESBatch[is.na(ESBatch)]=0
    } # end of if
    GS.WeightedBatch= as.matrix(datKMEBatch)^oddPower %*%  as.matrix(ESBatch)
    GS.Weighted=GS.Weighted+GS.WeightedBatch
  } # end of for (i in 1:nBlocks
  if (nMEs-nBlocks*blockSize>0 ) 
  {
    restindex=c((nBlocks*blockSize+1):nMEs)
    datMEBatch= datME[,restindex]
    datKMEBatch=as.matrix(signedKME(datExpr,datMEBatch, outputColumnName="MM."))
    ESBatch=   hubGeneSignificance(datKMEBatch ^oddPower,GS^oddPower)
    # the following omits the diagonal when datME=datExpr
    if (nGenes==nMEs) {diag(datKMEBatch[restindex,])=0
        # missing values will not be used 
        datKMEBatch[is.na(datKMEBatch)]=0
        ESBatch[is.na(ESBatch)]=0
    } # end of if (nGenes==nMEs) 
    GS.WeightedBatch= as.matrix(datKMEBatch)^oddPower %*% ESBatch
    GS.Weighted=GS.Weighted+GS.WeightedBatch
  } # end of if (nMEs-nBlocks*blockSize>0 )
  GS.Weighted=GS.Weighted/nMEs
  GS.Weighted[nAvailable< minimumSampleSize]=NA

  rankGS.Weighted=rank(-GS.Weighted, ties.method="first")
  rankGS=rank(-GS, ties.method="first")
  printFlush(paste("Proportion of agreement between GS.Weighted and GS:"))
  for (i in c(10,20,50,100,200,500,1000)) 
  {
    printFlush(paste("Top ", i, " list of genes: prop. of agreement = ", 
                signif(sum(rankGS.Weighted<=i & rankGS<=i,na.rm=T)/i,3)   ))
  } # end of for loop
  if (mean(abs(GS.Weighted),na.rm=T)>0) 
  {
    GS.Weighted=GS.Weighted/mean(abs(GS.Weighted),na.rm=T)*mean(abs(GS),na.rm=T)
  }
  if (addGS ) GS.Weighted=apply(data.frame(GS.Weighted, GS), 1,mean, na.rm=T)
  datout=data.frame(GS.Weighted, GS)

  datout
} # end of function

#--------------------------------------------------------------------------------------------------
#
# networkScreening
#
#--------------------------------------------------------------------------------------------------

networkScreening = function(
               y, datME, datExpr, 
               corFnc = "cor", corOptions = "use = 'p'",
               oddPower = 3,
               blockSize = 1000,
               minimumSampleSize = ..minNSamples,
               addMEy = TRUE, removeDiag = FALSE, 
               weightESy=0.5,
               getQValues = TRUE)
{
  oddPower=as.integer(oddPower)
  if (as.integer(oddPower/2)==oddPower/2 ) {oddPower=oddPower+1}
  nMEs=dim(as.matrix(datME))[[2]]
  nGenes=dim(as.matrix(datExpr))[[2]]
  # Here we add y as extra ME
  if (nGenes>nMEs & addMEy) {   datME=data.frame(y,datME)  }
  nMEs=dim(as.matrix(datME))[[2]]
  RawCor.Weighted=rep(0,nGenes)
  #Cor.Standard= as.numeric(cor(y,datExpr,use= "p") )
  corExpr = parse(text = paste("as.numeric( ", corFnc, "(y,datExpr ", prepComma(corOptions), "))")); 
  Cor.Standard= eval(corExpr)

  NoAvailable=apply(!is.na(datExpr), 2,sum)
  Cor.Standard[NoAvailable< minimumSampleSize]=NA
  if (nGenes==1) 
  {
    #RawCor.Weighted=as.numeric(cor(y,datExpr,use= "p") )
    corExpr = parse(text = paste("as.numeric(" , corFnc, "(y,datExpr ", prepComma(corOptions), "))"));
    RawCor.Weighted = eval(corExpr);
  }
  start = 1; i = 1; 
  while (start <= nMEs)
  {
    end = min(start + blockSize -1, nMEs);
    if (i>1 || end < nMEs) printFlush(paste("block number = ", i))
    index1=c(start:end)
    datMEBatch= datME[,index1]
    datKMEBatch=as.matrix(signedKME(datExpr,datMEBatch, outputColumnName="MM.", 
                                    corFnc = corFnc, corOptions = corOptions))
    # ES.CorBatch= as.vector(cor(  as.numeric(as.character(y))  ,datMEBatch, use="p"))
    corExpr = parse(text = paste("as.vector( ", corFnc, "(  as.numeric(as.character(y))  ,datMEBatch",
                                  prepComma(corOptions), "))" ));
    ES.CorBatch = eval(corExpr);

    #weightESy
    ES.CorBatch[ES.CorBatch>.999]= weightESy*1+ (1- weightESy)* 
                                    max(abs(ES.CorBatch[ES.CorBatch <.999 ]),na.rm=T)
    # the following omits the diagonal when datME=datExpr
    if (nGenes==nMEs & removeDiag) {diag(datKMEBatch[index1,])=0}
    if (nGenes==nMEs )
    {
      # missing values will not be used 
      datKMEBatch[is.na(datKMEBatch)]=0
      ES.CorBatch[is.na(ES.CorBatch)]=0
    } # end of if
    RawCor.WeightedBatch= as.matrix(datKMEBatch)^oddPower %*%  as.matrix(ES.CorBatch^oddPower)
    RawCor.Weighted=RawCor.Weighted+RawCor.WeightedBatch
    start = end + 1;
  } # end of while (start <= nMEs)
  RawCor.Weighted=RawCor.Weighted/nMEs
  RawCor.Weighted[NoAvailable< minimumSampleSize]=NA
  #to avoid dividing by zero we scale it as follows
  if (max(abs(RawCor.Weighted),na.rm=T)==1) RawCor.Weighted=RawCor.Weighted/1.0000001
  if (max(abs( Cor.Standard),na.rm=T)==1)  Cor.Standard=Cor.Standard/1.0000001
  RawZ.Weighted=sqrt(NoAvailable -2)*RawCor.Weighted/sqrt(1-RawCor.Weighted^2)
  Z.Standard= sqrt(NoAvailable -2)* Cor.Standard/sqrt(1-Cor.Standard^2)
  
  if (sum(abs(Z.Standard),na.rm=T) >0 ) 
  {
    Z.Weighted=RawZ.Weighted/sum(abs(RawZ.Weighted),na.rm=T)*sum(abs(Z.Standard),na.rm=T)
  } # end of if 
  h1=Z.Weighted/sqrt(NoAvailable-2)
  Cor.Weighted=h1/sqrt(1+h1^2)
  p.Weighted=as.numeric(2*(1-pt(abs(Z.Weighted),NoAvailable-2)))
  p.Standard=2*(1-pt(abs(Z.Standard),NoAvailable-2))

  if (getQValues)
  {
    # since the function qvalue cannot handle missing data, we set missing p-values to 1.
    p.Weighted2=p.Weighted
    p.Standard2=p.Standard
    p.Weighted2[is.na(p.Weighted)]=1
    p.Standard2[is.na(p.Standard)]=1
    
    q.Weighted=try(qvalue(p.Weighted2)$qvalues, silent = TRUE)
    q.Standard=try(qvalue(p.Standard2)$qvalues, silent = TRUE)
  
    if (class(q.Weighted)=="try-error") 
    {
      warning("Calculation of weighted q-values failed; the q-values will be returned as NAs.");
      q.Weighted=rep(NA, length(p.Weighted) )
    }
    if (class(q.Standard)=="try-error")
    {
      warning("Calculation of standard q-values failed; the q-values will be returned as NAs.");
      q.Standard=rep(NA, length(p.Standard) )
    }
  } else {
    q.Weighted=rep(NA, length(p.Weighted) )
    q.Standard=rep(NA, length(p.Standard) )
    if (getQValues)
      printFlush("networkScreening: Warning: package qvalue not found. q-values will not be calculated.");
  }
  rankCor.Weighted=rank(-abs(Cor.Weighted), ties.method="first")
  rankCor.Standard=rank(-abs(Cor.Standard), ties.method="first")
  printFlush(paste("Proportion of agreement between lists based on abs(Cor.Weighted) and abs(Cor.Standard):"))
  for (i in c(10,20,50,100,200,500,1000)) 
  {
    printFlush(paste("Top ", i, " list of genes: prop. agree = ", 
                signif(sum(rankCor.Weighted<=i & rankCor.Standard<=i,na.rm=T)/i,3)))
  } # end of for loop


  datout=data.frame(p.Weighted, q.Weighted, Cor.Weighted, Z.Weighted,
                    p.Standard, q.Standard, Cor.Standard, Z.Standard)
  names(datout) = sub("Cor", corFnc, names(datout), fixed = TRUE);
  datout
} # end of function


##############################################################################################
#
# Functions included from NetworkFunctions-PL-07.R
# Selected ones only
#
##############################################################################################


#---------------------------------------------------------------------------------------------------------
# labeledHeatmap.R
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#
# .reverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


.reverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

.reverseVector = function(Vector)
{
  ind = seq(from=length(Vector), to=1, by=-1);
  Vector[ind];
  #Vector
}
  
#--------------------------------------------------------------------------
#
# labeledHeatmap = function ( Matrix, xLabels, yLabels, ... ) { 
#
#--------------------------------------------------------------------------
# This function plots a heatmap of the specified matrix 
# and labels the x and y axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent 
# with the dimensions in.
# If colorLabels==TRUE, the labels are not printed and instead interpreted as colors --
#  -- a simple symbol with the appropriate color is printed instead of the label.
# The x,yLabels are expected to have the form "..color" as in "MEgrey" or "PCturquoise".
# xSymbol, ySymbols are additional markers that can be placed next to color labels

labeledHeatmap = function (Matrix, xLabels, yLabels = NULL, 
                           xSymbols = NULL, ySymbols = NULL, 
                           colorLabels = NULL, 
                           xColorLabels = FALSE, yColorLabels = FALSE,
                           checkColorsValid = TRUE,
                           invertColors = FALSE, 
                           setStdMargins = TRUE,
                           xLabelsPosition = "bottom",
                           xLabelsAngle = 45,
                           xLabelsAdj = 1,
                           xColorWidth = 0.05,
                           yColorWidth = 0.05,
                           colors = NULL, 
                           textMatrix = NULL, cex.text = NULL, cex.lab = NULL, 
                           cex.lab.x = cex.lab,
                           cex.lab.y = cex.lab,
                           colors.lab.x = 1,
                           colors.lab.y = 1,
                           plotLegend = TRUE, ... ) 
{
  if (!is.null(colorLabels)) {xColorLabels = colorLabels; yColorLabels = colorLabels; }
  
  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1]==dim(Matrix)[2])) 
    yLabels = xLabels; 

  if (checkColorsValid)
  {
    xValidColors = !is.na(match(substring(xLabels, 3), colors()));
    yValidColors = !is.na(match(substring(yLabels, 3), colors()));
  } else {
    xValidColors = rep(TRUE, length(xLabels));
    yValidColors = rep(TRUE, length(yLabels));
  }

  if (sum(xValidColors)>0) xColorLabInd = c(1:length(xLabels))[xValidColors]
  if (sum(!xValidColors)>0) xTextLabInd = c(1:length(xLabels))[!xValidColors]

  if (sum(yValidColors)>0) yColorLabInd = c(1:length(yLabels))[yValidColors]
  if (sum(!yValidColors)>0) yTextLabInd = c(1:length(yLabels))[!yValidColors]

  if (setStdMargins)
  {
    if (xColorLabels & yColorLabels)
    {
      par(mar=c(2,2,3,5)+0.2);
    } else {
      par(mar = c(7,7,3,5)+0.2);
    }
  }

  xLabPos = charmatch(xLabelsPosition, c("bottom", "top"));
  if (is.na(xLabPos))
    stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'");

  if (is.null(colors)) colors = heat.colors(30);
  if (invertColors) colors = .reverseVector(colors);
  labPos = .heatmapWithLegend(Matrix, signed = FALSE, colors = colors, cex.legend = cex.lab, 
                              plotLegend = plotLegend, ...)
  #if (plotLegend)
  #{
  #  image.plot(t(.reverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  #} else {
  #  image(z = t(.reverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  #}
  nxlabels = length(xLabels)
  plotbox = labPos$box;
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;

  xspacing = labPos$xMid[2] - labPos$xMid[1];
  yspacing = abs(labPos$yMid[2] - labPos$yMid[1]);

  nylabels = length(yLabels)
  offsetx = par("cxy")[1] / 3;
  offsety = par("cxy")[2] / 3;
  # Transform fractional widths into coordinate widths
  xColW = min(xmax - xmin, ymax - ymin) * xColorWidth;
  yColW = min(xmax - xmin, ymax - ymin) * yColorWidth;
  if (sum(!xValidColors)>0)
  {
    xLabYPos = ifelse(xLabPos==1, ymin - offsety, ymax + offsety)
    if (is.null(cex.lab)) cex.lab = 1;
    text(labPos$xMid[xTextLabInd] , xLabYPos, srt = xLabelsAngle, 
          adj = xLabelsAdj, labels = xLabels[xTextLabInd], xpd = TRUE, cex = cex.lab.x, col = colors.lab.x)
  }
  if (sum(xValidColors)>0)
  {
    baseY = ifelse(xLabPos==1, ymin-offsety-xColW, ymax + offsety + xColW);
    deltaY = ifelse(xLabPos==1, xColW, -xColW);
    rect(xleft = labPos$xMid[xColorLabInd] - xspacing/2, ybottom = baseY,
         xright = labPos$xMid[xColorLabInd] + xspacing/2, ytop = baseY + deltaY,
         density = -1,  col = substring(xLabels[xColorLabInd], 3), 
         border = substring(xLabels[xColorLabInd], 3), xpd = TRUE)
    if (!is.null(xSymbols))
      text ( labPos$xMid[xColorLabInd], baseY - sign(deltaY)* offsety, xSymbols[xColorLabInd], 
             adj = xLabelsAdj, 
             xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x, col = colors.lab.x);
  }
  if (sum(!yValidColors)>0)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(xmin - offsetx, labPos$yMid[yTextLabInd], srt = 0, 
         adj = c(1, 0.5), labels = yLabels[yTextLabInd], xpd = TRUE, cex = cex.lab.y, col = colors.lab.y )
  } 
  if (sum(yValidColors)>0)
  {
    rect(xleft = xmin- yColW - offsetx, ybottom = labPos$yMid[yColorLabInd] - yspacing/2,
         xright = xmin- offsetx, ytop = labPos$yMid[yColorLabInd] + yspacing/2, 
         density = -1,  col = substring(yLabels[yColorLabInd], 3), 
         border = substring(yLabels[yColorLabInd], 3), xpd = TRUE)
    if (!is.null(ySymbols))
      text (xmin- yColW - 2*offsetx, 
            labPos$yMid[yColorLabInd], ySymbols[yColorLabInd], 
            adj = c(1, 0.5), xpd = TRUE, cex = cex.lab.y, col = colors.lab.y);
  }

  if (!is.null(textMatrix))
  {
    if (is.null(cex.text)) cex.text = par("cex");
    if (!isTRUE(all.equal(dim(textMatrix), dim(Matrix))))
      stop("labeledHeatmap: textMatrix was given, but has dimensions incompatible with Matrix.");
    for (rw in 1:dim(Matrix)[1])
      for (cl in 1:dim(Matrix)[2])
      {
        text(labPos$xMid[cl], labPos$yMid[rw],
             as.character(textMatrix[rw,cl]), xpd = TRUE, cex = cex.text, adj = c(0.5, 0.5));
      }
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
}

#--------------------------------------------------------------------------
#
# labeledBarplot = function ( Matrix, labels, ... ) { 
#
#--------------------------------------------------------------------------
#
# Plots a barplot of the Matrix and writes the labels underneath such that they are readable.

labeledBarplot = function ( Matrix, labels, colorLabels = FALSE, colored = TRUE, 
                            setStdMargins = TRUE, stdErrors = NULL, cex.lab = NULL, 
                            xLabelsAngle = 45, ... ) 
{ 
  if (setStdMargins) par(mar=c(3,3,2,2)+0.2)

  if (colored)
  {
     colors = substring(labels, 3);
  } else {
     colors = rep("grey", times = ifelse(length(dim(Matrix))<2, length(Matrix), dim(Matrix)[[2]]));
  }

  ValidColors = !is.na(match(substring(labels, 3), colors()));
  
  if (sum(ValidColors)>0) ColorLabInd = c(1:length(labels))[ValidColors]
  if (sum(!ValidColors)>0) TextLabInd = c(1:length(labels))[!ValidColors]

  colors[!ValidColors] = "grey";
  
  mp = barplot(Matrix, col = colors, xaxt = "n", xlab="", yaxt="n", ...)

  if (length(dim(Matrix))==2) {
     means = apply(Matrix, 2, sum);
  } else {
     means = Matrix;  
  }

  if (!is.null(stdErrors)) addErrorBars(means, 1.96*stdErrors, two.side = TRUE);
  
  # axis(1, labels = FALSE)
  nlabels = length(labels)
  plotbox = par("usr");
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4];
  # print(paste("yrange:", yrange));
  if (nlabels>1)
  {
     spacing = (mp[length(mp)] - mp[1])/(nlabels-1);
  } else {
     spacing = (xmax-xmin);
  }
  yoffset = yrange/30
  xshift = spacing/2;
  xrange = spacing * nlabels;
  if (is.null(cex.lab)) cex.lab = 1;
  if (colorLabels)
  {
    #rect(xshift + ((1:nlabels)-1)*spacing - spacing/2.1, ymin - spacing/2.1 - spacing/8,
    #     xshift + ((1:nlabels)-1)*spacing + spacing/2.1, ymin - spacing/8,
    #     density = -1,  col = substring(labels, 3), border = substring(labels, 3), xpd = TRUE)
    if (sum(!ValidColors)>0)
    {
      text( mp[!ValidColors] , ymin - 0.02, srt = 45,
            adj = 1, labels = labels[TextLabInd], xpd = TRUE, cex = cex.lab, 
            srt = xLabelsAngle)
    }
    if (sum(ValidColors)>0)
    {
      rect(mp[ValidColors] - spacing/2.1, ymin - 2*spacing/2.1 * yrange/xrange - yoffset,
           mp[ValidColors] + spacing/2.1, ymin - yoffset,
           density = -1,  col = substring(labels[ValidColors], 3), 
           border = substring(labels[ValidColors], 3), xpd = TRUE)
    }
  } else {
    text(((1:nlabels)-1)*spacing +spacing/2 , ymin - 0.02*yrange, srt = 45, 
          adj = 1, labels = labels, xpd = TRUE, cex = cex.lab, srt = xLabelsAngle)
  }
  axis(2, labels = T)
}

#--------------------------------------------------------------------------
#
# sizeGrWindow
#
#--------------------------------------------------------------------------
# if the current device isn't of the required dimensions, close it and open a new one.

sizeGrWindow = function(width, height)
{
  din = par("din");
  if ( (din[1]!=width) | (din[2]!=height) )
  {
    dev.off();
    dev.new(width = width, height=height);
  }
}

#======================================================================================================
# GreenToRed.R
#======================================================================================================

greenBlackRed = function(n, gamma = 1)
{
  half = as.integer(n/2);
  red = c(rep(0, times = half), 0, seq(from=0, to=1, length.out = half)^(1/gamma));
  green = c(seq(from=1, to=0, length.out = half)^(1/gamma), rep(0, times = half+1));
  blue = rep(0, times = 2*half+1);
  col = rgb(red, green, blue, maxColorValue = 1);
  col;
}

greenWhiteRed = function(n, gamma = 1, warn = TRUE)
{
  if (warn) 
      warning(spaste("WGCNA::greenWhiteRed: this palette is not suitable for people\n",
                     "with green-red color blindness (the most common kind of color blindness).\n",
                     "Consider using the function blueWhiteRed instead."));
  half = as.integer(n/2);
  red = c(seq(from=0, to=1, length.out = half)^(1/gamma), rep(1, times = half+1));
  green = c(rep(1, times = half+1), seq(from=1, to=0, length.out = half)^(1/gamma));
  blue = c(seq(from=0, to=1, length.out = half)^(1/gamma), 1, 
          seq(from=1, to=0, length.out = half)^(1/gamma));
  col = rgb(red, green, blue, maxColorValue = 1);
  col;
}

redWhiteGreen = function(n, gamma = 1)
{
  half = as.integer(n/2);
  green = c(seq(from=0, to=1, length.out = half)^(1/gamma), rep(1, times = half+1));
  red = c(rep(1, times = half+1), seq(from=1, to=0, length.out = half)^(1/gamma));
  blue = c(seq(from=0, to=1, length.out = half)^(1/gamma), 1, 
               seq(from=1, to=0, length.out = half)^(1/gamma));
  col = rgb(red, green, blue, maxColorValue = 1);
  col;
}

#======================================================================================================
#
# Color pallettes that are more friendly to people with common color blindness
#
#======================================================================================================

blueWhiteRed = function(n, gamma = 1)
{
  blueEnd = c(0, 0.35, 0.80);
  redEnd = c(0.8, 0.2, 0);
  middle = c(1,1,1);

  half = as.integer(n/2);
  if (n%%2 == 0)
  {
    index1 = c(1:half);
    index2 = c(1:half)+half;
    frac1 = ((index1-1)/(half-1))^(1/gamma);
    frac2 = rev(frac1);
  } else {
    index1 = c(1:(half + 1))
    index2 = c(1:half) + half + 1
    frac1 = (c(0:half)/half)^(1/gamma);
    frac2 = rev((c(1:half)/half)^(1/gamma));
  }
  cols = matrix(0, n, 3);
  for (c in 1:3)
  {
    cols[ index1, c] = blueEnd[c] + (middle[c] - blueEnd[c]) * frac1;
    cols[ index2, c] = redEnd[c] + (middle[c] - redEnd[c]) * frac2;
  }

  rgb(cols[, 1], cols[, 2], cols[, 3], maxColorValue = 1);
}

#=========================================================================================================
#
# KeepCommonProbes
#
#-------------------------------------------------------------------------------------------
# Filters out probes that are not common to all datasets, and puts probes into the same order in each
# set. Works by creating dataframes of probe names and their indices and merging them all.

keepCommonProbes = function(multiExpr, orderBy = 1)
{
  size = checkSets(multiExpr);
  nSets = size$nSets;
  if (nSets<=0) stop("No expression data given!");

  Names = data.frame(Names = names(multiExpr[[orderBy]]$data));

  if (nSets>1) for (set in (1:nSets))
  {
    SetNames = data.frame(Names = names(multiExpr[[set]]$data), 
                          index = c(1:dim(multiExpr[[set]]$data)[2]));
    Names = merge(Names, SetNames, by.x = "Names", by.y = "Names", all = FALSE, sort = FALSE);
  }

  for (set in 1:nSets)
    multiExpr[[set]]$data = multiExpr[[set]]$data[, Names[, set+1]];

  multiExpr;
}
  
#--------------------------------------------------------------------------------------
#
# addTraitToPCs
#
#--------------------------------------------------------------------------------------

# Adds a trait vector to a set of eigenvectors.
# Caution: multiTraits is assumed to be a vector of lists with each list having an entry data which is 
# a nSamples x nTraits data frame with an appropriate column name, not a vector.

addTraitToMEs = function(multiME, multiTraits)
{
  nSets = length(multiTraits);
  setsize = checkSets(multiTraits);
  nTraits = setsize$nGenes;
  nSamples = setsize$nSamples;

  if (length(multiME)!=nSets)
    stop("Numbers of sets in multiME and multiTraits parameters differ - must be the same.");

  multiMETs = vector(mode="list", length=nSets);
  for (set in 1:nSets)
  {
    trait.subs = multiTraits[[set]]$data;
    multiMET = as.data.frame(cbind(multiME[[set]]$data, trait.subs));
    colnames(multiMET) = c(colnames(multiME[[set]]$data), colnames(trait.subs));
    if (!is.null(multiME[[set]]$AET))
    {
      AET = as.data.frame(cbind(multiME[[set]]$averageExpr, trait.subs));
      colnames(AET) = c(colnames(multiME[[set]]$averageExpr), colnames(trait.subs));
    }
    multiMETs[[set]] = list(data=multiMET);
  }
  multiMETs;
}


#--------------------------------------------------------------------------------------
#
# CorrelationPreservation
#
#--------------------------------------------------------------------------------------
#
# Given a set of multiME (or OrderedMEs), calculate the preservation values for each module in each pair
# of datasets and return them as a matrix

correlationPreservation = function(multiME, setLabels, excludeGrey = TRUE, greyLabel = "grey")
{
  nSets = length(multiME);
  if (nSets!=length(setLabels)) stop("The lengths of multiME and setLabels must equal.");
  if (nSets<=1) stop("Something is wrong with argument multiME: its length is 0 or 1");
  Names = names(multiME[[1]]$data);
  if (excludeGrey)
  {
      Use = substring(Names, 3)!=greyLabel;
  } else {
      Use = rep(TRUE, times = length(Names));
  }
  No.Mods = ncol(multiME[[1]]$data[, Use]); 
  CP = matrix(0, nrow = No.Mods, ncol = nSets*(nSets-1)/2);
  diag(CP) = 1;
  CPInd = 1;
  CPNames = NULL;
  for (i in 1:(nSets-1))
    for (j in (i+1):nSets)
    {
      corME1 = cor(multiME[[i]]$data[, Use], use="p");
      corME2 = cor(multiME[[j]]$data[, Use], use="p");
      d = 1-abs(tanh((corME1 - corME2) / (abs(corME1) + abs(corME2))^2));
      CP[ ,CPInd] = apply(d, 1, sum)-1;
      CPNames = c(CPNames, paste(setLabels[i], "::", setLabels[j], collapse = ""));
      CPInd = CPInd + 1;
    }
  CPx = as.data.frame(CP);
  names(CPx) = CPNames;
  rownames(CPx) = Names[Use];
  CPx;
}


#--------------------------------------------------------------------------------------
#
# setCorrelationPreservation
#
#--------------------------------------------------------------------------------------
#
# Given a set of multiME (or OrderedMEs), calculate the preservation values for each each pair
# of datasets and return them as a matrix.

setCorrelationPreservation = function(multiME, setLabels, excludeGrey = TRUE, greyLabel = "grey",
                                      method = "absolute")
{
  m = charmatch(method, c("absolute", "hyperbolic"));
  if (is.na(m))
  {
    stop("Unrecognized method given. Recognized methods are absolute, hyperbolic. ");
  }
  nSets = length(multiME);
  if (nSets!=length(setLabels)) stop("The lengths of multiME and setLabels must equal.");
  if (nSets<=1) stop("Something is wrong with argument multiME: its length is 0 or 1");
  Names = names(multiME[[1]]$data);
  if (excludeGrey)
  {
      Use = substring(Names, 3)!=greyLabel;
  } else {
      Use = rep(TRUE, times = length(Names));
  }
  No.Mods = ncol(multiME[[1]]$data[, Use]);
  SCP = matrix(0, nrow = nSets, ncol = nSets);
  diag(SCP) = 0;
  for (i in 1:(nSets-1))
    for (j in (i+1):nSets)
    {
      corME1 = cor(multiME[[i]]$data[, Use], use="p");
      corME2 = cor(multiME[[j]]$data[, Use], use="p");
      if (m==1) {
        d = 1 - abs(corME1 - corME2)/2;
      } else {
        d = 1-abs(tanh((corME1 - corME2) / (abs(corME1) + abs(corME2))^2));
      }
      SCP[i,j] = sum(d[upper.tri(d)])/sum(upper.tri(d));
      SCP[j,i] = SCP[i,j];
    }
  SCPx = as.data.frame(SCP);
  names(SCPx) = setLabels;
  rownames(SCPx) = setLabels;
  SCPx;
}

#---------------------------------------------------------------------------------------
#
# preservationNetworkDensity
#
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# 
# preservationNetworkConnectivity
#
#---------------------------------------------------------------------------------------

# This function returns connectivities of nodes in preservation networks

preservationNetworkConnectivity = function(
   multiExpr, 
   useSets = NULL, useGenes = NULL,
   corFnc = "cor", corOptions = "use='p'",
   networkType = "unsigned",
   power = 6,
   sampleLinks = NULL, nLinks = 5000,
   blockSize = 1000,
   setSeed = 12345,
   weightPower = 2,
   verbose = 2, indent = 0)

{
  spaces = indentSpaces(indent)

  size = checkSets(multiExpr);
  nGenes = size$nGenes;
  nSets = size$nSets;
  if (!is.null(useSets) || !is.null(useGenes))
  {
    if (is.null(useSets)) useSets = c(1:nSets)
    if (is.null(useGenes)) useGenes = c(1:nGenes)
    useExpr = vector(mode = "list", length = length(useSets));
    for (set in 1:length(useSets))
      useExpr[[set]] = list(data = multiExpr[[useSets[set]]]$data[, useGenes]);
    multiExpr = useExpr;
    rm(useExpr); collectGarbage();
  }
  size = checkSets(multiExpr);
  nGenes = size$nGenes;
  nSets = size$nSets;

  if (is.null(sampleLinks))
  {
    sampleLinks = (nGenes > nLinks);
  }

  if (sampleLinks) nLinks = min(nLinks, nGenes) else nLinks = nGenes;

  if (blockSize * nLinks > .largestBlockSize) blockSize = as.integer(.largestBlockSize/nLinks);

  intNetworkType = charmatch(networkType, .networkTypes);
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument. Recognized values are (unique abbreviations of)",
               paste(.networkTypes, collapse = ", ")));

  subtract = rep(1, nGenes);
  if (sampleLinks)
  {
    if (verbose > 0) 
      printFlush(paste(spaces, "preservationNetworkConnectivity: selecting sample pool of size",
                       nLinks, ".."))
    sd = apply(multiExpr[[1]]$data, 2, sd, na.rm = TRUE);
    order = order(-sd);
    saved = FALSE;
    if (exists(".Random.seed")) 
    {
      saved = TRUE;
      savedSeed = .Random.seed
      if (is.numeric(setSeed)) set.seed(setSeed);
    }
    samplePool = order[sample(x = nGenes, size = nLinks)]
    if (saved) .Random.seed <<- savedSeed;
    subtract[-samplePool] = 0;
  } 

  nPairComps = nSets * (nSets -1)/2;
  
  allPres = rep(NA, nGenes);
  allPresW = rep(NA, nGenes);
  allPresH = rep(NA, nGenes);
  allPresWH = rep(NA, nGenes);

  pairPres = matrix(NA, nGenes, nPairComps);
  pairPresW = matrix(NA, nGenes, nPairComps);
  pairPresH = matrix(NA, nGenes, nPairComps);
  pairPresWH = matrix(NA, nGenes, nPairComps);

  compNames = NULL;
  for (set1 in 1:(nSets-1))
    for (set2 in (set1+1):nSets)
      compNames = c(compNames, paste(set1, "vs", set2));

  dimnames(pairPres) = list(names(multiExpr[[1]]$data), compNames);
  dimnames(pairPresW) = list(names(multiExpr[[1]]$data), compNames);
  dimnames(pairPresH) = list(names(multiExpr[[1]]$data), compNames);
  dimnames(pairPresWH) = list(names(multiExpr[[1]]$data), compNames);

  if (verbose>0) 
  {
     pind = initProgInd(trailStr = " done");
  }

  nBlocks = as.integer((nGenes-1)/blockSize);
  SetRestrConn = NULL;
  start = 1;
  if (sampleLinks)
  {
    corEval = parse(text = paste(corFnc, 
                       "(multiExpr[[set]]$data[, samplePool], multiExpr[[set]]$data[, blockIndex] ", 
                       prepComma(corOptions), ")"))
  } else {
    corEval = parse(text = paste(corFnc, 
                       "(multiExpr[[set]]$data, multiExpr[[set]]$data[, blockIndex] ", 
                        prepComma(corOptions), ")"))
  }

  while (start <= nGenes)
  {
    end = start + blockSize-1;
    if (end>nGenes) end = nGenes;
    blockIndex = c(start:end);
    nBlockGenes = end-start+1;
    blockAdj = array(0, dim = c(nSets, nLinks, nBlockGenes));
    #if (verbose>1) printFlush(paste(spaces, "..working on genes", start, "through", end, "of", nGenes))
    for (set in 1:nSets)
    {
      c = eval(corEval);
      if (intNetworkType==1)
      { c = abs(c);
      } else if (intNetworkType==2)
      { c = (1+c)/2;
      } else if (intNetworkType==3)
      { c[c < 0] = 0;
      } else stop("Internal error: intNetworkType has wrong value:", intNetworkType, ". Sorry!");
      adj_mat = as.matrix(c^power);
      if (sum(is.na(adj_mat)) > 0)
        stop("NA values present in adjacency - this function cannot handle them yet. Sorry!");
      adj_mat[is.na(adj_mat)] = 0;
      blockAdj[set, , ] = adj_mat
    }
    min = matrix(0, nLinks, nBlockGenes)
    which = matrix(0, nLinks, nBlockGenes)
    res = .C("minWhichMin", as.double(blockAdj), as.integer(nSets), as.integer(nLinks * nBlockGenes),
                    min = as.double(min), as.double(which), DUP = FALSE)
    min[, ] = res$min;
    max = matrix(0, nLinks, nBlockGenes);
    res = .C("minWhichMin", as.double(-blockAdj), as.integer(nSets), as.integer(nLinks * nBlockGenes),
                    min = as.double(min), as.double(which), DUP = FALSE)
    max[, ] = -res$min;
    rm(res);
    diff = max - min;
    allPres[blockIndex] = (apply(1-diff, 2, sum) - subtract[blockIndex])/(nLinks - subtract[blockIndex]);
    weight = ((max + min)/2)^weightPower
    allPresW[blockIndex] = (apply((1-diff) * weight, 2, sum) - subtract[blockIndex])/
                              (apply(weight, 2, sum) - subtract[blockIndex]);
    hyp = 1-tanh(diff/(max+min)^2);
    allPresH[blockIndex] = (apply(hyp, 2, sum) - subtract[blockIndex])/(nLinks - subtract[blockIndex]);
    allPresWH[blockIndex] = (apply(hyp * weight, 2, sum) - subtract[blockIndex])/
                              (apply(weight, 2, sum) - subtract[blockIndex]);

    compNames = NULL;
    compInd = 1;
    for (set1 in 1:(nSets-1))
      for (set2 in (set1+1):nSets)
      {
        diff = abs(blockAdj[set1, , ] - blockAdj[set2, , ]) 
        compNames = c(compNames, paste(set1, "vs", set2));
        pairPres[blockIndex, compInd] = (apply(1-diff, 2, sum) - subtract[blockIndex]) /
                                        (nLinks - subtract[blockIndex]);
        weight = ((blockAdj[set1, , ] + blockAdj[set2, , ])/2)^weightPower
        pairPresW[blockIndex, compInd] = (apply((1-diff) * weight, 2, sum) - subtract[blockIndex]) /
                                        (apply(weight, 2, sum) - subtract[blockIndex]);
        hyp = 1-tanh(diff/(blockAdj[set1, , ] + blockAdj[set2, , ])^2)
        pairPresH[blockIndex, compInd] = (apply(hyp, 2, sum) - subtract[blockIndex]) /
                                        (nLinks - subtract[blockIndex]);
        pairPresWH[blockIndex, compInd] = (apply(hyp * weight, 2, sum) - subtract[blockIndex]) /
                                        (apply(weight, 2, sum) - subtract[blockIndex]);
        compInd = compInd + 1;
      }

    start = end+1;
    if (verbose>0) pind = updateProgInd(end/nGenes, pind);
    collectGarbage();
  }
  if (verbose>0) printFlush(" ");
  list(pairwise = pairPres, complete = allPres, pairwiseWeighted = pairPresW,
       completeWeighted = allPresW, pairwiseHyperbolic = pairPresH, completeHyperbolic = allPresH,
       pairwiseWeightedHyperbolic = pairPresWH, completeWeightedHyperbolic = allPresWH)
}

#--------------------------------------------------------------------------------------
#
# plotEigengeneNetworks
#
#--------------------------------------------------------------------------------------
# Plots a matrix plot of the ME(T)s. On the diagonal the heatmaps show correlation of MEs in the
# particular subset; off-diagonal are differences in the correlation matrix. 
# setLabels is a vector of titles for the diagonal diagrams; the off-diagonal will have no title
# for now.

plotEigengeneNetworks = function(
                      multiME, setLabels,
                      letterSubPlots = FALSE, Letters = NULL, 
                      excludeGrey = TRUE, greyLabel = "grey", 
                      plotDendrograms = TRUE,
                      plotHeatmaps = TRUE,
                      setMargins = TRUE, 
                      marDendro = NULL, marHeatmap = NULL, 
                      colorLabels = TRUE, signed = TRUE,
                      heatmapColors = NULL,
                      plotAdjacency = TRUE, 
                      printAdjacency = FALSE, cex.adjacency = 0.9,
                      coloredBarplot = TRUE, barplotMeans = TRUE, barplotErrors = FALSE,
                      plotPreservation = "standard", 
                      zlimPreservation = c(0,1),
                      printPreservation = FALSE, cex.preservation = 0.9, 
                      ...)
{
  # invertColors = FALSE;
  size = checkSets(multiME, checkStructure = TRUE);
  if (!size$structureOK)
  {
    #printFlush(paste(
    #  "plotEigengeneNetworks: Given multiME does not appear to be a multi-set structure.\n",
    #  "Will attempt to convert it into a multi-set structure containing 1 set."));
    multiME = fixDataStructure(multiME);
  }

  if (is.null(Letters)) Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  if (is.null(heatmapColors)) 
    if (signed)
    {
      heatmapColors = greenWhiteRed(50);
    } else {
      heatmapColors = heat.colors(30);
    }
  nSets = length(multiME);
  cex = par("cex");
  mar = par("mar");
  nPlotCols = nSets;
  nPlotRows = as.numeric(plotDendrograms) + nSets * as.numeric(plotHeatmaps);
  if (nPlotRows==0)
    stop("Nothing to plot: neither dendrograms not heatmaps requested.")
  par(mfrow = c(nPlotRows, nPlotCols));
  par(cex = cex);
  if (excludeGrey) for (set in 1:nSets)
    multiME[[set]]$data = 
        multiME[[set]]$data[ , substring(names(multiME[[set]]$data),3)!=greyLabel]

  plotPresTypes = c("standard", "hyperbolic", "both")
  ipp = pmatch(plotPreservation, plotPresTypes);
  if (is.na(ipp))
    stop(paste("Invalid 'plotPreservation'. Available choices are", 
               paste(plotPresTypes, sep = ", ")));
  
  letter.ind = 1;
  if (plotDendrograms) for (set in 1:nSets)
  {
    #par(cex = StandardCex/1.4);
    par(mar = marDendro);
    labels = names(multiME[[set]]$data);
    uselabels = labels[substring(labels,3)!=greyLabel];
    corME = cor(multiME[[set]]$data[substring(labels,3)!=greyLabel,
                                 substring(labels,3)!=greyLabel], use="p");
    disME = as.dist(1-corME);
    clust = flashClust(disME, method = "average");
    if (letterSubPlots) {
      main = paste(substring(Letters, letter.ind, letter.ind), ". ", setLabels[set], sep="");
    } else {
      main = setLabels[set];
    }
    #validColors = is.na(match(uselabels, colors()));
    #plotLabels = ifelse(validColors, substring(uselabels[validColors], 3), uselabels[!validColors]);
    plotLabels = uselabels;
    plot(clust, main = main, sub="", xlab="", 
         labels = plotLabels, ylab="", ylim=c(0,1));
    letter.ind = letter.ind + 1;
  }

  if (plotHeatmaps) for (i.row in (1:nSets)) for (i.col in (1:nSets))
  {
    letter.ind = i.row * nSets + i.col;
    if (letterSubPlots) 
    {
       #letter = paste("(", substring(Letters, first = letter.ind, last = letter.ind), ")", sep = "");
       letter = paste( substring(Letters, first = letter.ind, last = letter.ind), ".  ", sep = "");
    } else {
       letter = NULL;
    }
    par(cex = cex);
    if (setMargins) {
      if (is.null(marHeatmap))
      {
        if (colorLabels) {
          par(mar = c(1,2,3,4)+0.2);
        } else {
          par(mar = c(6,7,3,5)+0.2);
        }
      } else {
        par(mar = marHeatmap);
      }
    }
    nModules = dim(multiME[[i.col]]$data)[2]
    textMat = NULL;
    if (i.row==i.col)
    {
      corME = cor(multiME[[i.col]]$data, use="p") 
      pME = corPvalueFisher(corME, nrow(multiME[[i.col]]$data));
      if (printAdjacency)
      {
         textMat = paste(signif(corME, 2), "\n", signif(pME, 1));
         dim(textMat) = dim(corME)
      } 
      if (signed)
      {
        if (plotAdjacency) {
         if (printAdjacency) 
         {
            textMat = paste(signif((1+corME)/2, 2), "\n", signif(pME, 1));
            dim(textMat) = dim(corME)
         } 
         labeledHeatmap((1+corME)/2, names(multiME[[i.col]]$data), names(multiME[[i.col]]$data),
                               main=paste(letter, setLabels[[i.col]]), invertColors=FALSE,
                               zlim=c(0,1.0),
                               colorLabels = colorLabels, colors = heatmapColors, 
                               setStdMargins = FALSE, 
                               textMatrix = textMat, cex.text = cex.adjacency, ...);
        } else {
         labeledHeatmap(corME, names(multiME[[i.col]]$data), names(multiME[[i.col]]$data),
                               main=paste(letter, setLabels[[i.col]]), invertColors=FALSE,
                               zlim=c(-1,1.0),
                               colorLabels = colorLabels, colors = heatmapColors, setStdMargins = FALSE, 
                               textMatrix = textMat, cex.text = cex.adjacency, ...);
        }
      } else {
         labeledHeatmap(abs(corME), names(multiME[[i.col]]$data), names(multiME[[i.col]]$data),
                               main=paste(letter, setLabels[[i.col]]), invertColors=FALSE,
                               zlim=c(0,1.0),
                               colorLabels = colorLabels, colors = heatmapColors, 
                               setStdMargins = FALSE, 
                               textMatrix = textMat, cex.text = cex.adjacency, ...);
      }
    } else
    {
      corME1 = cor(multiME[[i.col]]$data, use="p");
      corME2 = cor(multiME[[i.row]]$data, use="p");
      cor.dif = (corME1 - corME2)/2;
      d = tanh((corME1 - corME2) / (abs(corME1) + abs(corME2))^2);
      # d = abs(corME1 - corME2) / (abs(corME1) + abs(corME2));
      if (ipp==1 | ipp==3) 
      {
         dispd = cor.dif;
         main = paste(letter, "Preservation");
         if (ipp==3) {
            dispd[upper.tri(d)] = d[upper.tri(d)];
            main=paste(letter, "Hyperbolic preservation (UT)\nStandard preservation (LT)")
         }
      } else {
         dispd = d;
         main = paste(letter, "Hyperbolic preservation");
      }
      if (i.row>i.col)
      {
        if (signed)
        {
          half = as.integer(length(heatmapColors)/2);
          range = c(half:length(heatmapColors)); 
          halfColors = heatmapColors[range];
        } else {
          halfColors = heatmapColors;
        }
        if (printPreservation) {
          printMtx = matrix(paste(".", as.integer((1-abs(dispd))*100), sep = ""), 
                             nrow = nrow(dispd), ncol = ncol(dispd));
          printMtx[printMtx==".100"] = "1";
        } else { 
          printMtx = NULL; 
        }
        if ( (sum( (1-abs(dispd))<zlimPreservation[1]) || ((1-abs(dispd))>zlimPreservation[2])) >0)
          warning("plotEigengeneNetworks: Correlation preservation data out of zlim range.");
        labeledHeatmap(1-abs(dispd), names(multiME[[i.col]]$data), names(multiME[[i.col]]$data), 
                       main = main, invertColors=FALSE,
                       colorLabels = colorLabels, zlim = zlimPreservation, colors = halfColors,
                       setStdMargins = FALSE, 
                       textMatrix = printMtx, cex.text = cex.preservation, ...);
      } else {
        if (ipp==2) {
           dp = 1-abs(d);
           method = "Hyperbolic:";
        } else {
           dp = 1-abs(cor.dif); 
           method = "Preservation:";
        }
        diag(dp) = 0;
        if (barplotMeans) {
          sum_dp = mean(dp[upper.tri(dp)]);
          means = apply(dp, 2, sum)/(ncol(dp)-1);
          if (barplotErrors) {
             errors = sqrt( (apply(dp^2, 2, sum)/(ncol(dp)-1) - means^2)/(ncol(dp)-2));
          } else {
             errors = NULL; 
          }
          labeledBarplot(means, names(multiME[[i.col]]$data), 
                         main=paste(letter, "D=", signif(sum_dp,2)), 
                         ylim=c(0,1),
                         colorLabels = colorLabels, colored = coloredBarplot,
                         setStdMargins = FALSE, stdErrors = errors, ... )
        } else {
          sum_dp = sum(dp[upper.tri(dp)]);
          labeledBarplot(dp, names(multiME[[i.col]]$data),
                         main=paste(letter, method, "sum = ", signif(sum_dp,3)), 
                         ylim=c(0,dim(dp)[[1]]),
                         colorLabels = colorLabels, colored = coloredBarplot, 
                         setStdMargins = FALSE, ... )
        }
      }
    }
  }
}

#====================================================================================================
#
# numbers2colors: convert a vector of numbers to colors
#
#====================================================================================================

# Turn a numerical variable into a color indicator. x can be a matrix or a vector.
# For discrete variables, consider also labels2colors.

numbers2colors = function(x, 
                     signed = NULL, 
                     centered = signed,
                     lim = NULL, 
                     commonLim = FALSE,
                     colors = if (signed) blueWhiteRed(100) else blueWhiteRed(100)[51:100],
                     naColor = "grey")
{
  x = as.matrix(x);
  if (!is.numeric(x))
    stop("'x' must be numeric. For a factor, please use as.numeric(x) in the call.");
  if (is.null(signed))
  {
     if (any(x<0, na.rm = TRUE) & any(x>0, na.rm = TRUE))
     {
       signed = TRUE;
     } else
       signed = FALSE;
  }
  if (is.null(centered)) centered = signed;

  if (is.null(lim))
  {
    if (signed & centered)
    {
      max = apply(abs(x), 2, max, na.rm = TRUE);
      lim = as.matrix(cbind(-max, max));
    } else {
      lim = as.matrix(cbind(apply(x, 2, min, na.rm = TRUE),  apply(x, 2, max, na.rm = TRUE)));
    }
    if (commonLim) 
      lim = c(min(lim[, 1], na.rm = TRUE), max(lim[, 2], na.rm = TRUE));
  }
  if (is.null(dim(lim)))
  {
    if (length(lim)!=2)
      stop("'lim' must be a vector of length 2 or a matrix with 2 columns.");
    if (!is.numeric(lim))
      stop("'lim' must be numeric");
    if (sum(is.finite(lim))!=2) stop("'lim' must be finite.");
    lim = t(as.matrix(lim));
  } else {
    if (ncol(x)!=nrow(lim))
      stop("Incompatible numbers of columns in 'x' and rows in 'lim'.")
    if (!is.numeric(lim))
      stop("'lim' must be numeric");
    if (sum(is.finite(lim))!=length(lim)) stop("'lim' must be finite.");
  }

  xMin = matrix(lim[,1], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  xMax = matrix(lim[,2], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  if (sum(xMin==xMax)>0)
    warning("(some columns in) 'x' are constant. Their color will be the color of NA.");

  xx = x;
  xx[is.na(xx)] = ((xMin+xMax)[is.na(xx)])/2;
  if (sum(x < xMin, na.rm = TRUE) > 0)
  {
    warning("Some values of 'x' are below given minimum and will be truncated to the minimum.");
    x[xx<xMin] = xMin[xx<xMin];
  }

  if (sum(x > xMax, na.rm = TRUE) > 0)
  {
    warning("Some values of 'x' are above given maximum and will be truncated to the maximum.");
    x[xx>xMax] = xMax[xx>xMax];
  }

  mmEq = xMin==xMax;

  nColors = length(colors);

  xCol = array(naColor, dim = dim(x));

  xInd = (x - xMin)/(xMax-xMin);
  xInd[xInd==1] = 1-0.5/nColors;
  xCol[!mmEq] = colors[as.integer(xInd[!mmEq] * nColors) + 1];
  xCol[is.na(xCol)] = naColor;
  
  xCol;
}

#====================================================================================================
#
#  Rand index calculation
#
#====================================================================================================

# this function is used for computing the Rand index below...
#
.choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {out1[i] <- 0}
    else {out1[i] <- choose(n[i], k)}}
  out1	
}


# the following function computes the Rand index between 2 clusterings
# assesses how similar two clusterings are
randIndex <- function(tab, adjust=TRUE) 
{
  a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
  m <- nrow(tab);
  n <- ncol(tab);
  for (i in 1:m) {
    c<-0
    for(j in 1:n) {
      a <- a+.choosenew(tab[i,j],2)
      nj <- sum(tab[,j])
      c <- c+.choosenew(nj,2)
    }
    ni <- sum(tab[i,])
    b <- b+.choosenew(ni,2)
    nn <- nn+ni
  }
  if(adjust==T) {
    d <- .choosenew(nn,2)
    adrand <- (a-(b*c)/d)/(0.5*(b+c)-(b*c)/d)
    adrand
  } else {
    b <- b-a
    c <- c-a
    d <- .choosenew(nn,2)-a-b-c
    rand <- (a+d)/(a+b+c+d)
    rand
  }
}


#============================================================================================
#
# Check expression data: mark genes and samples with too many missing entries
#
#============================================================================================

goodGenes = function(datExpr, useSamples = NULL, useGenes = NULL,
                     minFraction = 1/2, minNSamples = ..minNSamples, minNGenes = ..minNGenes,
                     verbose = 1, indent = 0)
{
  if (is.null(useGenes)) useGenes = rep(TRUE, ncol(datExpr));
  if (is.null(useSamples)) useSamples = rep(TRUE, nrow(datExpr));

  if (length(useGenes)!= ncol(datExpr))
    stop("Length of nGenes is not compatible with number of columns in datExpr.");
  if (length(useSamples)!= nrow(datExpr))
    stop("Length of nSamples is not compatible with number of rows in datExpr.");

  nSamples = sum(useSamples);
  nGenes = sum(useGenes);
  nPresent = apply(!is.na(datExpr[useSamples, useGenes]), 2, sum)
  gg = useGenes;
  gg[useGenes][nPresent<minNSamples] = FALSE;
  var = apply(datExpr[useSamples, gg], 2, sd, na.rm = TRUE);
  nNAsGenes = apply(is.na(datExpr[useSamples, gg]), 2, sum);
  gg[gg] = (nNAsGenes < (1-minFraction) * nSamples & var>0 & (nSamples-nNAsGenes >= minNSamples));
  if (sum(gg) < minNGenes)
    stop("Too few genes with valid expression levels in the required number of samples.");

  if (verbose>0 & (nGenes - sum(gg) > 0))
    printFlush(paste("  ..Excluding", nGenes - sum(gg),
                     "genes from the calculation due to too many missing samples or zero variance."));

  gg;
}

goodSamples = function(datExpr, useSamples = NULL, useGenes = NULL,
                     minFraction = 1/2, minNSamples = ..minNSamples, minNGenes = ..minNGenes,
                     verbose = 1, indent = 0)
{
  if (is.null(useGenes)) useGenes = rep(TRUE, ncol(datExpr));
  if (is.null(useSamples)) useSamples = rep(TRUE, nrow(datExpr));

  if (length(useGenes)!= ncol(datExpr))
    stop("Length of nGenes is not compatible with number of columns in datExpr.");
  if (length(useSamples)!= nrow(datExpr))
    stop("Length of nSamples is not compatible with number of rows in datExpr.");

  nSamples = sum(useSamples);
  nGenes = sum(useGenes);
  nNAsSamples = apply(is.na(datExpr[useSamples, useGenes]), 1, sum);
  goodSamples = useSamples;
  goodSamples[useSamples] = ((nNAsSamples < (1-minFraction)*nGenes) & 
                             (nGenes - nNAsSamples >= minNGenes));
  if (sum(goodSamples) < minNSamples)
    stop("Too few samples with valid expression levels for the required number of genes.");

  if (verbose>0 & (nSamples - sum(goodSamples)>0))
    printFlush(paste("  ..Excluding", nSamples - sum(goodSamples),
                     "samples from the calculation due to too many missing genes."));

  goodSamples;
}

goodGenesMS = function(multiExpr, useSamples = NULL, useGenes = NULL,
                       minFraction = 1/2, minNSamples = ..minNSamples, minNGenes = ..minNGenes,
                       verbose = 1, indent = 0)
{
  dataSize = checkSets(multiExpr);
  nSets = dataSize$nSets;
  if (is.null(useGenes)) useGenes = rep(TRUE, dataSize$nGenes);
  if (is.null(useSamples)) 
  {
    useSamples = list();
    for (set in 1:nSets) useSamples[[set]] = rep(TRUE, dataSize$nSamples[set]);
  }

  if (length(useGenes)!= dataSize$nGenes)
    stop("Length of nGenes is not compatible with number of genes in multiExpr.");
  if (length(useSamples)!= nSets)
    stop("Length of nSamples is not compatible with number of sets in multiExpr.");

  for (set in 1:nSets) if (length(useSamples[[set]])!=dataSize$nSamples[set])
    stop(paste("Number of samples in useSamples[[", set, "]] incompatible\n   ",
               "with number of samples in the corresponding set of multiExpr."))

  nSamples = sapply(useSamples, sum);
  nGenes = sum(useGenes);

  goodGenes = useGenes;
  for (set in 1:nSets)
  {
    if (sum(goodGenes)==0) break;
    if (sum(useSamples[[set]])==0) next;
    nPresent = apply(!is.na(multiExpr[[set]]$data[useSamples[[set]], goodGenes]), 2, sum)
    goodGenes[goodGenes] = (nPresent >= minNGenes)
    var = apply(multiExpr[[set]]$data[useSamples[[set]], goodGenes], 2, sd, na.rm = TRUE);
    nNAsGenes = apply(is.na(multiExpr[[set]]$data[useSamples[[set]], goodGenes]), 2, sum);
    goodGenes[goodGenes][nNAsGenes > (1-minFraction)*nSamples[set] | var==0 | 
                           (nSamples[set]-nNAsGenes < minNSamples)] = FALSE;
  }
  if (sum(goodGenes) < minNGenes)
    stop("Too few genes with valid expression levels in the required number of samples in all sets.");

  if (verbose>0 & (nGenes - sum(goodGenes) > 0))
    printFlush(paste("  ..Excluding", nGenes - sum(goodGenes),
                     "genes from the calculation due to too many missing samples or zero variance."));
  goodGenes;
}

goodSamplesMS = function(multiExpr, useSamples = NULL, useGenes = NULL,
                         minFraction = 1/2, minNSamples = ..minNSamples, minNGenes = ..minNGenes,
                         verbose = 1, indent = 0)
{
  dataSize = checkSets(multiExpr);
  nSets = dataSize$nSets;
  if (is.null(useGenes)) useGenes = rep(TRUE, dataSize$nGenes);
  if (is.null(useSamples))
  {
    useSamples = list();
    for (set in 1:nSets) useSamples[[set]] = rep(TRUE, dataSize$nSamples[set]);
  }

  if (length(useGenes)!= dataSize$nGenes)
    stop("Length of nGenes is not compatible with number of genes in multiExpr.");
  if (length(useSamples)!= dataSize$nSets)
    stop("Length of nSamples is not compatible with number of sets in multiExpr.");

  for (set in 1:nSets) if (length(useSamples[[set]])!=dataSize$nSamples[set])
    stop(paste("Number of samples in useSamples[[", set, "]] incompatible\n   ",
               "with number of samples in the corresponding set of multiExpr."))

  nSamples = sapply(useSamples, sum);
  nGenes = sum(useGenes);

  goodSamples = useSamples;
  for (set in 1:nSets)
  {
    if (sum(useGenes)==0) break;
    if (sum(goodSamples[[set]])==0) next;
    nNAsSamples = apply(is.na(multiExpr[[set]]$data[useSamples[[set]], useGenes]), 1, sum);
    goodSamples[[set]][useSamples[[set]]] = 
          ((nNAsSamples < (1-minFraction) * nGenes) & (nGenes - nNAsSamples >= minNGenes));
    if (sum(goodSamples[[set]]) < minNSamples)
      stop("Too few samples with valid expression levels for the required number of genes in set", set);
    if (verbose>0 & (nSamples[set] - sum(goodSamples[[set]])>0))
      printFlush(paste("  ..Set", set,": Excluding", nSamples[set] - sum(goodSamples[[set]]),
                       "samples from the calculation due to too many missing genes."));
  }
  goodSamples;
}

goodSamplesGenes = function(datExpr, minFraction = 1/2, minNSamples = ..minNSamples, 
                            minNGenes = ..minNGenes, verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent)
  goodGenes = NULL;
  goodSamples = NULL;
  nBadGenes = 0;
  nBadSamples = 0;
  changed = TRUE;
  iter = 1;
  if (verbose>0)
      printFlush(paste(spaces, "Flagging genes and samples with too many missing values..."));
  while (changed)
  {
    if (verbose>0)
      printFlush(paste(spaces, " ..step", iter));
    goodGenes = goodGenes(datExpr, goodSamples, goodGenes,
                            minFraction = minFraction, minNSamples = minNSamples,
                            minNGenes = minNGenes, verbose = verbose - 1, indent = indent + 1);
    goodSamples = goodSamples(datExpr, goodSamples, goodGenes,
                            minFraction = minFraction, minNSamples = minNSamples,
                            minNGenes = minNGenes, verbose = verbose - 1, indent = indent + 1);
    changed = ( (sum(!goodGenes)>nBadGenes) | (sum(!goodSamples)>nBadSamples) )
    nBadGenes = sum(!goodGenes);
    nBadSamples = sum(!goodSamples);
    iter = iter + 1;
  }
  allOK = (sum(c(nBadGenes, nBadSamples)) == 0)
  list(goodGenes = goodGenes, goodSamples = goodSamples, allOK = allOK);
}

goodSamplesGenesMS = function(multiExpr, minFraction = 1/2, minNSamples = ..minNSamples, 
                              minNGenes = ..minNGenes, verbose = 2, indent = 0)
{
  spaces = indentSpaces(indent)
  size = checkSets(multiExpr)
  nSets = size$nSets;
  goodGenes = NULL;
  goodSamples = NULL;
  nBadGenes = 0;
  nBadSamples = rep(0, nSets);
  changed = TRUE;
  iter = 1;
  if (verbose>0)
      printFlush(paste(spaces, "Flagging genes and samples with too many missing values..."));
  while (changed)
  {
    if (verbose>0)
      printFlush(paste(spaces, " ..step", iter));
    goodGenes = goodGenesMS(multiExpr, goodSamples, goodGenes,
                            minFraction = minFraction, minNSamples = minNSamples,
                            minNGenes = minNGenes, verbose = verbose - 1, indent = indent + 1);
    goodSamples = goodSamplesMS(multiExpr, goodSamples, goodGenes,
                            minFraction = minFraction, minNSamples = minNSamples,
                            minNGenes = minNGenes, verbose = verbose - 1, indent = indent + 1);
    changed = FALSE;
    for (set in 1:nSets)
      changed = ( changed | (sum(!goodGenes)>nBadGenes) | (sum(!goodSamples[[set]])>nBadSamples[set]) )
    nBadGenes = sum(!goodGenes);
    for (set in 1:nSets) nBadSamples[set] = sum(!goodSamples[[set]]);
    iter = iter + 1;
    if (verbose > 2) 
       printFlush(paste(spaces, "   ..bad gene count: ", nBadGenes, 
                        ", bad sample counts: ", paste(nBadSamples, collapse = ", "), sep=""));
  }
  allOK = (sum(c(nBadGenes, nBadSamples)) == 0)
  list(goodGenes = goodGenes, goodSamples = goodSamples, allOK = allOK);
}

#============================================================================================
#
# modified heatmap plot: allow specifying the hang parameter for both side and top dendrograms
#
#============================================================================================
.heatmap = function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
    distfun = dist, hclustfun = flashClust, reorderfun = function(d, 
        w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
        "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
    margins = c(1.2, 1.2), ColSideColors, RowSideColors, cexRow = 0.2 + 
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
    verbose = getOption("verbose"), setLayout = TRUE, hang = 0.04, ...) 
{
    scale <- if(symm && missing(scale)) "none" else match.arg(scale)
    if(length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1L]
    nc <- di[2L]
    if(nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")
    if(!is.numeric(margins) || length(margins) != 2L)
        stop("'margins' must be a numeric vector of length 2")

    doRdend <- !identical(Rowv,NA)
    doCdend <- !identical(Colv,NA)
    if(!doRdend && identical(Colv, "Rowv")) doCdend <- FALSE
    ## by default order by row/col means
    if(is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
    if(is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)

    ## get the dendrograms and reordering indices
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            if (class(hcr)=='hclust')
            {
              hcr$height = hcr$height-min(hcr$height) + hang * (max(hcr$height)-min(hcr$height));
            }
            ddr <- as.dendrogram(hcr, hang = hang)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) x else t(x)))
            if (class(hcr)=='hclust')
            {
              hcc$height = hcc$height-min(hcc$height) + hang * (max(hcc$height)-min(hcc$height));
            }
            ddc <- as.dendrogram(hcc, hang = hang)
            if (!is.logical(Colv) || Colv) ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc

    ## reorder x
    x <- x[rowInd, colInd];

    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }

    ## Calculate the plot layout
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.5 else 0, 4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, "; lmat=\n")
        print(lmat)
    }
    if (!symm || scale != "none") x <- t(x)
    op <- par(no.readonly = TRUE)
    if (revC) {
        iy <- nc:1
        ddr <- rev(ddr)
        rowInd.colors = rev(rowInd)
        x <- x[, iy]
    } else iy <- 1:nr
    #on.exit(par(op))
    # print(paste("main:", main));
    if (setLayout) layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd.colors], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, 
          xlab = "", ylab = "", ...)
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
    if (!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    if (!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none" )
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1.8 else 0, margins[2]))
    if (doCdend) plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none" )
    else if (!is.null(main)) frame()
    if (!is.null(main)) title(main, cex.main = 1.2 * op[["cex.main"]])
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}


#===================================================================================================
# The vectorize functions turns a matrix or data frame into a vector. If the matrix is not #symmetric the
# number of entries of the vector equals the number of rows times the #number of columns of the matrix.
# But if the matrix is symmetrical then it only uses the #entries in the upper triangular matrix.
#If the option diag =TRUE, it also includes the diagonal elements of the symmetric #matrix. By default it
# excludes the diagonal elements of a symmetric matrix.

vectorizeMatrix=function(M, diag=FALSE) 
{
  if ( is.null(dim(M)) )  stop("The input of the vectorize function is not a matrix or data frame.")
  if ( length(dim(M))!=2 )  stop("The input of the vectorize function is not a matrix or data frame.")
  # now we check whether the matrix is symmetrical
  if (dim(M)[[1]]==dim(M)[[2]])
  {
      M=as.matrix(M)
      Mtranspose=t(M)
      abs.difference=max( abs(M-Mtranspose),na.rm=T)
      if (abs.difference<10^(-14) ) 
      {
          out=M[upper.tri(M,diag)]  
      }
      else
          out=as.vector(M);
  } else
      out=as.vector(M)
  out
} # end

#========================================================================================================

scaleFreeFitIndex=function(k,nBreaks=10, removeFirst = FALSE)
{
        discretized.k = cut(k, nBreaks)
        dk = tapply(k, discretized.k, mean)
        p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
        breaks1 = seq(from = min(k), to = max(k), 
            length = nBreaks + 1)
        hist1 = hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
        dk2 = hist1$mids
        dk = ifelse(is.na(dk), dk2, dk)
        dk = ifelse(dk == 0, dk2, dk)
        p.dk = ifelse(is.na(p.dk), 0, p.dk)
        log.dk = as.vector(log10(dk))
        if (removeFirst) {
            p.dk = p.dk[-1]
            log.dk = log.dk[-1]
        }
       log.p.dk= as.numeric(log10(p.dk + 1e-09))
        lm1 = lm(log.p.dk ~ log.dk)
        lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
   datout=data.frame(Rsquared.SFT=summary(lm1)$r.squared,
                     slope.SFT=summary(lm1)$coefficients[2, 1], 
                     truncatedExponentialAdjRsquared= summary(lm2)$adj.r.squared)
   datout
} # end of function scaleFreeFitIndex

#========================================================================================================

standardScreeningCensoredTime= function (
   time,
   event,
   datExpr,
   percentiles = seq(from = 0.1, to = 0.9, by = 0.2),
   dichotomizationResults = FALSE,
   qValues = TRUE,
   fastCalculation = TRUE)
{
datExpr=data.frame(datExpr)
    no.Columns = dim(as.matrix(datExpr))[[2]]
    m = dim(as.matrix(datExpr))[[1]]
    if (length(time) != m) 
        stop("The length of the time variable does not equal the number of rows of datExpr.\nConsider transposing datExpr.")
    if (length(event) != m) 
        stop("The length of the event variable does not equal the number of rows of datExpr.\nConsider transposing datExpr.")
    if (fastCalculation) {
        fittemp = summary(coxph(Surv(time, event) ~ 1, na.action = na.exclude))
        CumHazard = predict(fittemp, type = "expected")
        martingale1 = event - CumHazard
        deviance0 = ifelse(event == 0, 2 * CumHazard, -2 * log(CumHazard) + 
            2 * CumHazard - 2)
        devianceresidual = sign(martingale1) * sqrt(deviance0)
        corDeviance = as.numeric(cor(devianceresidual, datExpr, 
            use = "p"))
        no.nonMissing = sum(!is.na(time))
        pvalueDeviance = corPvalueFisher(cor = corDeviance, nSamples = no.nonMissing)
       qvalueDeviance=rep(NA, length(pvalueDeviance) )
                   rest1= ! is.na( pvalueDeviance) 
          qvalueDeviance [rest1] = qvalue(pvalueDeviance [rest1])$qvalues

        datout = data.frame(ID = dimnames(datExpr)[[2]], pvalueDeviance, 
            qvalueDeviance, corDeviance)
    }
    if (!fastCalculation) {
        pvalueWald = rep(NA, no.Columns)
        HazardRatio = rep(NA, no.Columns)
        CI.UpperLimitHR = rep(NA, no.Columns)
        CI.LowerLimitHR = rep(NA, no.Columns)
        C.index = rep(NA, no.Columns)
        pvalueLogrank = rep(NA, no.Columns)
        pValuesDichotomized = data.frame(matrix(NA, nrow = no.Columns, 
            ncol = length(percentiles)))
        names(pValuesDichotomized) = paste("pValueDichotPercentile", 
            as.character(percentiles), sep = "")
        fittemp = summary(coxph(Surv(time, event) ~ 1, na.action = na.exclude))
        CumHazard = predict(fittemp, type = "expected")
        martingale1 = event - CumHazard
        deviance0 = ifelse(event == 0, 2 * CumHazard, -2 * log(CumHazard) + 
            2 * CumHazard - 2)
        devianceresidual = sign(martingale1) * sqrt(deviance0)
        corDeviance = as.numeric(cor(devianceresidual, datExpr, 
            use = "p"))
        no.nonMissing = sum(!is.na(time))
        pvalueDeviance = corPvalueFisher(cor = corDeviance, nSamples = no.nonMissing)
      

for (i in 1:no.Columns) {
            Column = as.numeric(as.matrix(datExpr[, i]))
            var1 = var(Column, na.rm = T)
            if (var1 == 0 | is.na(var1)) {
                pvalueWald[i] = NA
                pvalueLogrank[i] = NA
                HazardRatio[i] = NA
                CI.UpperLimitHR[i] = NA
                CI.LowerLimitHR[i] = NA
                C.index[i] = NA
                 }  # end of              if (var1 == 0 | is.na(var1))
            if (var1 != 0 & !is.na(var1)) {
                cox1 = summary(coxph(Surv(time, event) ~ Column, 
                  na.action = na.exclude))
                pvalueWald[i] = cox1$coef[5]
                pvalueLogrank[i] = cox1$sctest[[3]]
                HazardRatio[i] = exp(cox1$coef[1])
                CI.UpperLimitHR[i] = exp(cox1$coef[1] + 1.96 * 
                  cox1$coef[3])
                CI.LowerLimitHR[i] = exp(cox1$coef[1] - 1.96 * 
                  cox1$coef[3])
                C.index[i] = rcorr.cens(Column, Surv(time, event), 
                  outx = T)[[1]]
            } # end of   if (var1 != 0 & !is.na(var1)) 


            if (dichotomizationResults) {
                quantilesE = as.numeric(quantile(Column, prob = percentiles))
                for (j in 1:length(quantilesE)) {
                  ColumnDichot = I(Column > quantilesE[j])
                  var1 = var(ColumnDichot, na.rm = T)
                  if (var1 == 0 | is.na(var1)) {
                    pValuesDichotomized[i, j] = NA
                  } # end of if
                  if (var1 != 0 & !is.na(var1)) {
                    coxh = summary(coxph(Surv(time, event) ~ 
                      ColumnDichot, na.action = na.exclude))
                    pValuesDichotomized[i, j] = coxh$coef[5]
                  } # end of if
                } # end of for (j
                MinimumDichotPvalue = apply(pValuesDichotomized, 
                  1, min, na.rm = T)
               } # end of if (dichotomizationResults)
            


if (!qValues) {
                datout = data.frame(ID = dimnames(datExpr)[[2]], 
                  pvalueWald, pvalueLogrank, pvalueDeviance, 
                  corDeviance, HazardRatio, CI.LowerLimitHR, 
                  CI.UpperLimitHR, C.index)
            }      # end of      if (!qValues) {

        } # end of for (i in 1:no.Columns) 


  if (qValues) {
       qvalueWald=rep(NA, length(pvalueWald) )
                   rest1= ! is.na( pvalueWald) 
          qvalueWald [rest1] = qvalue(pvalueWald[rest1])$qvalues

       qvalueLogrank=rep(NA, length(pvalueLogrank) )
                   rest1= ! is.na( pvalueLogrank) 
          qvalueLogrank [rest1] = qvalue(pvalueLogrank[rest1])$qvalues

       qvalueDeviance=rep(NA, length(pvalueDeviance) )
                   rest1= ! is.na( pvalueDeviance) 
          qvalueDeviance [rest1] = qvalue(pvalueDeviance[rest1])$qvalues

                datout = data.frame(ID = dimnames(datExpr)[[2]], 
                  pvalueWald, qvalueWald, pvalueLogrank, qvalueLogrank, 
                  pvalueDeviance,       qvalueDeviance , corDeviance, HazardRatio, CI.LowerLimitHR, 
                  CI.UpperLimitHR, C.index)
            } # end of  if (qValues)


        if (dichotomizationResults) {
            datout = data.frame(datout, MinimumDichotPvalue, 
                pValuesDichotomized)
        }
    }
    datout
} # end of function standardScreeningCensoredTime


#================================================================================
#
# standardScreeningNumericTrait
#
#================================================================================

standardScreeningNumericTrait= function (datExpr, yNumeric, corFnc = cor, 
                                         corOptions = list(use = 'p'),
                                         alternative = c("two.sided", "less", "greater"),
                                         qValues = TRUE, 
                                         areaUnderROC = TRUE) 
{ 
  datExpr=as.matrix(datExpr)
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  if (length(yNumeric) != nSamples)
      stop("the length of the sample trait y does not equal the number of rows of datExpr")
  corPearson = rep(NA, nGenes)
  pvalueStudent = rep(NA, nGenes);
  AreaUnderROC = rep(NA, nGenes);
  nPresent = Z = rep(NA, nGenes);
       
  corFnc = match.fun(corFnc);
  corOptions$y = yNumeric;
  corOptions$x = as.matrix(datExpr);
  cp = do.call(corFnc, corOptions);
  corPearson = as.numeric(cp);

  finMat = !is.na(datExpr)
  np = t(finMat) %*% (!is.na(as.matrix(yNumeric)))

  nPresent = as.numeric(np)

  ia = match.arg(alternative)
  T = sqrt(np - 2) * corPearson/sqrt(1 - corPearson^2)
  if (ia == "two.sided") {
       p = 2 * pt(abs(T), np - 2, lower.tail = FALSE)
  }
  else if (ia == "less") {
      p = pt(T, np - 2, lower.tail = TRUE)
  }
  else if (ia == "greater") {
      p = pt(T, np - 2, lower.tail = FALSE)
  }
  pvalueStudent = as.numeric(p);

  Z = 0.5 * log( (1+corPearson)/(1-corPearson) ) * sqrt(nPresent -2 );

  if (areaUnderROC) for (i in 1:dim(datExpr)[[2]]) 
  {
    AreaUnderROC[i] = rcorr.cens(datExpr[, i], yNumeric, outx = T)[[1]]
  }

  q.Student=rep(NA, length(pvalueStudent) )
  rest1= ! is.na(pvalueStudent) 
  if (qValues)
  {
    x = try({ q.Student[rest1] = qvalue(pvalueStudent[rest1])$qvalues }, silent = TRUE)
    if (inherits(x, "try-error"))
      printFlush(paste("Warning in standardScreeningNumericTrait: function qvalue returned an error.\n",
                       "The returned qvalues will be invalid. The qvalue error: ", x, "\n"));
  }
  if (is.null(colnames(datExpr)))
  {
     ID = spaste("Variable.", 1:ncol(datExpr));
  } else
     ID = colnames(datExpr);

  output = data.frame(ID = ID, cor = corPearson,
                      Z = Z,
                      pvalueStudent = pvalueStudent);
  if (qValues) output$qvalueStudent = q.Student;
  if (areaUnderROC) output$AreaUnderROC = AreaUnderROC;

  output$nPresentSamples = nPresent;

  output
}




#================================================================================
#
# spaste
#
#================================================================================

spaste = function(...) { paste(..., sep = "") }

#================================================================================
#
# metaZfunction
#
#================================================================================

metaZfunction=function(datZ, columnweights=NULL  )
{
  if ( ! is.null(columnweights) )  {datZ=   t(t(datZ)* columnweights)   } 
  datZpresent= !is.na(datZ)+0.0
  if ( ! is.null(columnweights) )  {datZpresent=   t(t(datZpresent)* columnweights)   } 
  sumZ=as.numeric(apply(datZ,1,sum,na.rm=T))
  variance= as.numeric(apply(datZpresent^2 ,1,sum))
  sumZ/sqrt(variance)
}

#================================================================================
#
# rankPvalue
#
#================================================================================

rankPvalue=function(datS, columnweights = NULL, na.last = "keep", ties.method = "average", 
    calculateQvalue = TRUE, pValueMethod = "all") 
{
    no.rows = dim(datS)[[1]]
    no.cols = dim(datS)[[2]]
    if (!is.null(columnweights) & no.cols != length(columnweights)) 
        stop("The number of components of the vector columnweights is unequal to the number of columns of datS. Hint: consider transposing datS. ")

if (!is.null(columnweights) ) {
if ( min(columnweights,na.rm=T)<0 )  stop("At least one component of columnweights is negative, which makes no sense. The entries should be positive numbers")
if ( sum(is.na(columnweights))>0 )  stop("At least one component of columnweights is missing, which makes no sense. The entries should be positive numbers")
if ( sum( columnweights)!= 1 ) {
 # warning("The entries of columnweights do not sum to 1. Therefore, they will divided by the sum. Then the resulting weights sum to 1.");
columnweights= columnweights/sum( columnweights)
}
}

    if (pValueMethod != "scale") {
              percentilerank1 = function(x) {
            R1 = rank(x, ties.method = ties.method, na.last = na.last)
            (R1-.5)/max(R1, na.rm = T)
        }
 
        datrankslow = apply(datS, 2, percentilerank1)
        if (!is.null(columnweights)) {
            datrankslow = t(t(datrankslow) * columnweights)
        }
        datSpresent = !is.na(datS) + 0
        if (!is.null(columnweights)) {
            datSpresent = t(t(datSpresent) * columnweights)
        }
        expectedsum = apply(datSpresent, 1, sum, na.rm = T) * 
            0.5
        varsum = apply(datSpresent^2, 1, sum, na.rm = T) * 1/12
        observed.sumPercentileslow = as.numeric(apply(datrankslow, 
            1, sum, na.rm = T))
        Zstatisticlow = (observed.sumPercentileslow - expectedsum)/sqrt(varsum)
        datrankshigh = apply(-datS, 2, percentilerank1)
        if (!is.null(columnweights)) {
            datrankshigh = t(t(datrankshigh) * columnweights)
        }
        observed.sumPercentileshigh = as.numeric(apply(datrankshigh, 
            1, sum, na.rm = T))
        Zstatistichigh = (observed.sumPercentileshigh - expectedsum)/sqrt(varsum)
        pValueLow = pnorm((Zstatisticlow))
        pValueHigh = pnorm((Zstatistichigh))
        pValueExtreme = pmin(pValueLow, pValueHigh)
        datoutrank = data.frame(pValueExtreme, pValueLow, pValueHigh)
        if (calculateQvalue) {
            qValueLow = rep(NA, dim(datS)[[1]])
            qValueHigh = rep(NA, dim(datS)[[1]])
            qValueExtreme = rep(NA, dim(datS)[[1]])
            rest1 = !is.na(pValueLow)
            qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
            rest1 = !is.na(pValueHigh)
            qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
            rest1 = !is.na(pValueExtreme)
            qValueExtreme = pmin(qValueLow, qValueHigh)
            datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
            datoutrank = data.frame(datoutrank, datq)
            names(datoutrank) = paste(names(datoutrank), "Rank", 
                sep = "")
        }
    }
    if (pValueMethod != "rank") {
        datSpresent = !is.na(datS) + 0
        scaled.datS = scale(datS)
        if (!is.null(columnweights)) {
            scaled.datS = t(t(scaled.datS) * columnweights)
            datSpresent = t(t(datSpresent) * columnweights)
        }
        expected.value = rep(0, no.rows)
        varsum = apply(datSpresent^2, 1, sum, na.rm = T) * 1
        observed.sumScaleddatS = as.numeric(apply(scaled.datS, 
            1, sum, na.rm = T))
        Zstatisticlow = (observed.sumScaleddatS - expected.value)/sqrt(varsum)
        scaled.minusdatS = scale(-datS)
        if (!is.null(columnweights)) {
            scaled.minusdatS = t(t(scaled.minusdatS) * columnweights)
        }
        observed.sumScaledminusdatS = as.numeric(apply(scaled.minusdatS, 
            1, sum, na.rm = T))
        Zstatistichigh = (observed.sumScaledminusdatS - expected.value)/sqrt(varsum)
        pValueLow = pnorm((Zstatisticlow))
        pValueHigh = pnorm((Zstatistichigh))
        pValueExtreme = 2 * pnorm(-abs(Zstatisticlow))
        datoutscale = data.frame(pValueExtreme, pValueLow, pValueHigh)
        if (calculateQvalue) {
            qValueLow = rep(NA, dim(datS)[[1]])
            qValueHigh = rep(NA, dim(datS)[[1]])
            qValueExtreme = rep(NA, dim(datS)[[1]])
            rest1 = !is.na(pValueLow)
            qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
            rest1 = !is.na(pValueHigh)
            qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
            rest1 = !is.na(pValueExtreme)
            qValueExtreme[rest1] = qvalue(pValueExtreme[rest1])$qvalues
            datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
            datoutscale = data.frame(datoutscale, datq)
        }
        names(datoutscale) = paste(names(datoutscale), "Scale", 
            sep = "")
    }
    if (pValueMethod == "rank") {
        datout = datoutrank
    }
    if (pValueMethod == "scale") {
        datout = datoutscale
    }
    if (pValueMethod != "rank" & pValueMethod != "scale") 
        datout = data.frame(datoutrank, datoutscale)
    datout
} # End of function

#========================================================================================================
#
# utility function: add a comma to string if the string is non-empty
#
#========================================================================================================

prepComma = function(s)
{
  if (s=="") return (s);
  paste(",", s);
}


#========================================================================================================
#
# "restricted" q-value calculation
#
#========================================================================================================

qvalue.restricted = function(p, trapErrors = TRUE, ...)
{
  fin = is.finite(p);
  qx = try(qvalue(p[fin], ...)$qvalues, silent = TRUE);
  q = rep(NA, length(p));
  if (inherits(qx, "try-error"))
  {
    if (!trapErrors) stop(qx);
  } else 
    q[fin] = qx;
  q;
}


#========================================================================================================
#
# consensusKME
#
#========================================================================================================

consensusKME = function(multiExpr, moduleLabels, multiEigengenes = NULL, consensusQuantile = 0,
                        signed = TRUE,
                        useModules = NULL,
                        metaAnalysisWeights = NULL, 
                        corAndPvalueFnc = corAndPvalue, corOptions = list(),
                        corComponent = "cor", getQvalues = FALSE,
                        useRankPvalue = TRUE,
                        rankPvalueOptions = list(calculateQvalue = getQvalues, pValueMethod = "scale"),
                        setNames = NULL, excludeGrey = TRUE,
                        greyLabel = ifelse(is.numeric(moduleLabels), 0, "grey"))
{
  corAndPvalueFnc = match.fun(corAndPvalueFnc);

  size = checkSets(multiExpr);
  nSets = size$nSets;
  nGenes = size$nGenes;
  nSamples = size$nSamples;

  if (!is.null(metaAnalysisWeights))
     if (length(metaAnalysisWeights)!=nSets)
       stop("Length of 'metaAnalysisWeights' must equal number of input sets.");

  if (!is.null(useModules))
  {
    if (greyLabel %in% useModules) 
      stop(paste("Grey module (or module 0) cannot be used with 'useModules'.\n",
                 "   Use 'excludeGrey = FALSE' to obtain results for the grey module as well. "));
    keep = moduleLabels %in% useModules;
    if (sum(keep)==0)
      stop("Incorrectly specified 'useModules': no such module(s).");
    moduleLabels [ !keep ] = greyLabel;
  }

  if (is.null(multiEigengenes))
    multiEigengenes = multiSetMEs(multiExpr, universalColors = moduleLabels, verbose = 0, 
                                  excludeGrey = excludeGrey, grey = greyLabel);

  modLevels = substring(colnames(multiEigengenes[[1]]$data), 3);
  nModules = length(modLevels);

  kME = p = Z = nObs = array(NA, dim = c(nGenes, nModules, nSets));

  corOptions$alternative = c("two.sided", "greater")[signed+1];
  
  haveZs = FALSE;
  for (set in 1:nSets)
  {
    corOptions$x = multiExpr[[set]]$data;
    corOptions$y = multiEigengenes[[set]]$data;
    cp = do.call(corAndPvalueFnc, args = corOptions);
    corComp = grep(corComponent, names(cp));
    pComp = match("p", names(cp));
    if (is.na(pComp)) pComp = match("p.value", names(cp));
    if (is.na(pComp)) stop("Function `corAndPvalueFnc' did not return a p-value.");
    kME[, , set] = cp[[corComp]]
    p[, , set] = cp[[pComp]];
    if (!is.null(cp$Z)) { Z[, , set] = cp$Z; haveZs = TRUE}
    if (!is.null(cp$nObs)) 
    {
       nObs[, , set] = cp$nObs;
    } else
       nObs[, , set] = t(is.na(multiExpr[[set]]$data)) %*% (!is.na(multiEigengenes[[set]]$data));
  }

  if (getQvalues)
  {
    q = apply(p, c(2:3), qvalue.restricted);
  } else q = NULL;

  # kME.average = rowMeans(kME, dims = 2); <-- not neccessary since weighted average also contains it

  powers = c(0, 0.5, 1);
  nPowers = length(powers)
  nWeights = nPowers + !is.null(metaAnalysisWeights)
  weightNames = c("equalWeights", "RootDoFWeights", "DoFWeights", "userWeights") [1:nWeights];
  kME.weightedAverage = array(NA, dim = c(nGenes, nWeights, nModules));
  for (m in 1:nWeights)
  {
    if (m<=nPowers) {
      weights = nObs^powers[m]
    } else
      weights = array( rep(metaAnalysisWeights, rep(nGenes*nModules, nSets)),
                             dim = c(nGenes, nModules, nSets));
    kME.weightedAverage[, m, ] = rowSums( kME * weights, na.rm = TRUE, dims = 2) / 
                                    rowSums(weights, dims = 2, na.rm = TRUE)
  }

  dim(kME.weightedAverage) = c(nGenes * nWeights, nModules);

  if (any(is.na(kME)))
  {
     kME.consensus.1 = apply(kME, c(1,2), quantile, prob = consensusQuantile, na.rm = TRUE);
     kME.consensus.2 = apply(kME, c(1,2), quantile, prob = 1-consensusQuantile, na.rm = TRUE);
     kME.median = apply(kME, c(1,2), median, na.rm = TRUE);
  } else {
    kME.consensus.1 = matrix( 
                      colQuantileC(t(matrix(kME, nGenes * nModules, nSets)), p = consensusQuantile),
                      nGenes, nModules);
    kME.consensus.2 = matrix( 
                      colQuantileC(t(matrix(kME, nGenes * nModules, nSets)), p = 1-consensusQuantile),
                      nGenes, nModules);
    kME.median = matrix(colQuantileC(t(matrix(kME, nGenes * nModules, nSets)), p = 0.5),
                        nGenes, nModules);
  }
  kME.consensus = ifelse(kME.median > 0, kME.consensus.1, kME.consensus.2);

  kME.consensus[ kME.consensus * kME.median < 0 ] = 0;

  # Prepare identifiers for the variables (genes)
  if (is.null(colnames(multiExpr[[1]]$data)))
  {
     ID = spaste("Variable.", 1:nGenes);
  } else
     ID = colnames(multiExpr[[1]]$data);

  # Get meta-Z, -p, -q values
  if (haveZs)
  {
    Z.kME.meta = p.kME.meta = array(0, dim = c(nGenes, nWeights, nModules))
    if (getQvalues) q.kME.meta = array(0, dim = c(nGenes, nWeights, nModules));
    for (m in 1:nWeights)
    {
      if (m<=nPowers) {
        weights = nObs^powers[m]
      } else
        weights = array( rep(metaAnalysisWeights, rep(nGenes*nModules, nSets)), 
                             dim = c(nGenes, nModules, nSets));

      Z1 = rowSums( Z * weights, na.rm = TRUE, dims = 2) / sqrt(rowSums(weights^2, na.rm = TRUE, dims = 2))
      if (signed)
      {
         p1 = pnorm(Z1, lower.tail = FALSE);
      } else
         p1 = 2*pnorm(abs(Z1), lower.tail = FALSE);
      Z.kME.meta[, m, ] = Z1;
      p.kME.meta[, m, ] = p1;
      if (getQvalues)
      {
        q1 = apply(p1, 2, qvalue.restricted);
        q.kME.meta[, m, ] = q1;
      }
    }
    dim(Z.kME.meta) = dim(p.kME.meta) = c(nGenes* nWeights, nModules);
    if (getQvalues) 
    {
        dim(q.kME.meta) = c(nGenes * nWeights, nModules);
    } else 
        q.kME.meta = NULL;
  } else {
    Z.kME.meta = p.kME.meta = q.kME.meta = NULL;
  }

  # Call rankPvalue

  if (useRankPvalue)
  {
    for (mod in 1:nModules) for (m in 1:nWeights)
    {
      if (m<=nPowers) {
        weights = nObs[, mod, ]^powers[m]
      } else
        weights = matrix( metaAnalysisWeights, nGenes, nSets, byrow = TRUE);
      # rankPvalue requires a vector of weights... so compress the weights to a vector.
      # Output a warning if the compression loses information.
      nDifferent = apply(weights, 2, function(x) {length(unique(x)) });
      if (any(nDifferent)>1)
        printFlush(paste("Warning in consensusKME: rankPvalue requires compressed weights.\n",
                         "Some weights may not be entirely accurate."));
      cw = colMeans(weights, na.rm = TRUE);
      rankPvalueOptions$columnweights = cw / sum(cw);

      rankPvalueOptions$datS = kME[, mod, ];
      rp1 = do.call(rankPvalue, rankPvalueOptions);
      colnames(rp1) = spaste(colnames(rp1), ".ME", modLevels[mod], ".", weightNames[m]);
      if (mod==1 && m==1) {
        rp = rp1;
      } else 
        rp = cbind(rp, rp1);
    }
  }

  # Format the output... this will entail some rearranging of the individual set results.
  if (is.null(setNames))
     setNames = names(multiExpr);

  if (is.null(setNames))
     setNames = spaste("Set_", c(1:nSets));

  if (!haveZs) Z = NULL;

  keep = c(TRUE, TRUE, getQvalues, haveZs);
  varNames = c("kME", "p.kME", "q.kME", "Z.kME")[keep];
  nVars = sum(keep);
  combined = array(c (kME, p, q, Z), dim = c(nGenes, nModules, nSets, nVars));

  recast = matrix( c(cast(melt(combined), X1~X4~X3~X2)), nGenes, nSets * nModules * nVars);
  colnames(recast) = spaste( rep(varNames, nSets * nModules),
                             rep( modLevels, rep(nVars * nSets, nModules)), 
                             ".", rep( rep(setNames, rep(nVars, nSets)), nModules));
  combinedMeta.0 = rbind(
             kME.consensus,
             kME.weightedAverage,
             Z.kME.meta,
             p.kME.meta,
             q.kME.meta);

  combinedMeta = matrix(combinedMeta.0, nGenes, 
                            (1 + nWeights + (2*haveZs + haveZs*getQvalues)*nWeights) * nModules);
  metaNames = c("consensus.kME", 
                spaste("weightedAverage.", weightNames, ".kME"), 
                spaste("meta.Z.", weightNames, ".kME"), 
                spaste("meta.p.", weightNames, ".kME"),
                spaste("meta.q.", weightNames, ".kME")
                )[ c(rep(TRUE, nWeights+1), rep(haveZs, nWeights), rep(haveZs, nWeights), 
                               rep(haveZs && getQvalues, nWeights))];
  nMetaVars = length(metaNames);
  colnames(combinedMeta) = spaste (rep(metaNames, nModules), 
                                   rep(modLevels, rep(nMetaVars, nModules)));

  if (useRankPvalue) {
     out = data.frame(ID = ID, combinedMeta, rp, recast);
  } else 
     out = data.frame(ID = ID, combinedMeta, recast);

  out
}

#======================================================================================================
#
# Meta-analysis
#
#======================================================================================================

.isBinary = function(multiTrait)
{
  bin = TRUE;
  for (set in 1:length(multiTrait))
    if (length(sort(unique(multiTrait[[set]]$data))) > 2) bin = FALSE;

  bin;
}

metaAnalysis = function(multiExpr, multiTrait, 
                        binary = NULL,
                        #consensusQuantile = 0,
                        metaAnalysisWeights = NULL,
                        corFnc = cor, corOptions = list(use = 'p'),
                        getQvalues = FALSE,
                        getAreaUnderROC = FALSE,
                        useRankPvalue = TRUE,
                        rankPvalueOptions = list(),
                        setNames = NULL, 
                        kruskalTest = FALSE, var.equal = FALSE, 
                        metaKruskal = kruskalTest,
                        na.action = "na.exclude")
{

  size = checkSets(multiExpr);
  nSets = size$nSets;

  for (set in 1:nSets)
    multiTrait[[set]] $ data = as.matrix(multiTrait[[set]] $ data);

  tSize = checkSets(multiTrait);
  if (tSize$nGenes!=1)
     stop("This function only works for a single trait. ");

  if (size$nSets!=tSize$nSets)
     stop("The number of sets in 'multiExpr' and 'multiTrait' must be the same.");

  if (!all.equal(size$nSamples, tSize$nSamples))
     stop("Numbers of samples in each set of 'multiExpr' and 'multiTrait' must be the same.");

  #if (!is.finite(consensusQuantile) || consensusQuantile < 0 || consensusQuantile > 1)
  #   stop("'consensusQuantile' must be between 0 and 1.");

  if (is.null(setNames))
     setNames = names(multiExpr);

  if (is.null(setNames))
     setNames = spaste("Set_", c(1:nSets));

  if (metaKruskal && !kruskalTest) 
     stop("Kruskal statistic meta-analysis requires kruskal test. Use kruskalTest=TRUE.");

  if (is.null(binary)) binary = .isBinary(multiTrait);

  if (!is.null(metaAnalysisWeights))
  {
    if (length(metaAnalysisWeights)!=nSets)
      stop("Length of 'metaAnalysisWeights' must equal the number of sets in 'multiExpr'.")
    if (any (!is.finite(metaAnalysisWeights)) || any(metaAnalysisWeights < 0))
      stop("All weights in 'metaAnalysisWeights' must be positive.");
  }

  setResults = list();

  for (set in 1:size$nSets)
  {
    if (binary)
    {
      setResults[[set]] = standardScreeningBinaryTrait(multiExpr[[set]]$data,
                            as.vector(multiTrait[[set]]$data), kruskalTest = kruskalTest, 
                            qValues = getQvalues, var.equal = var.equal, na.action = na.action,
                            corFnc = corFnc, corOptions = corOptions);
      trafo = TRUE;
      if (metaKruskal) 
      {
        metaStat = "stat.Kruskal.signed";
        metaP = "pvaluekruskal";
      } else {
        metaStat = "t.Student";
        metaP = "pvalueStudent"
      }
    } else {
      setResults[[set]] = standardScreeningNumericTrait(multiExpr[[set]]$data,
                            as.vector(multiTrait[[set]]$data), qValues = getQvalues, 
                            corFnc = corFnc, corOptions = corOptions, 
                            areaUnderROC = getAreaUnderROC);
      metaStat = "Z";
      trafo = FALSE;
    }
  }

  comb = NULL;
  for (set in 1:nSets)
  {
    if (set==1) 
    {
      comb =  setResults[[set]] [, -1];
      ID = setResults[[set]] [, 1];
      colNames= colnames(comb);
      nColumns = ncol(comb);
      colnames(comb) = spaste("X", c(1:nColumns));
    } else {
      xx = setResults[[set]][, -1];
      colnames(xx) = spaste("X", c(1:nColumns));
      comb = rbind(comb, xx);
    }
  }

  # Re-arrange comb:

  comb = matrix(as.matrix(as.data.frame(comb)), size$nGenes, nColumns * nSets);

  colnames(comb) = spaste( rep( colNames, rep(nSets, nColumns)), ".", rep(setNames, nColumns));

  # Find the columns from which to do meta-analysis
  statCols = grep(spaste("^", metaStat), colnames(comb));
  if (length(statCols)==0) stop("Internal error: no columns for meta-analysis found. Sorry!");
  setStats = comb[, statCols];

  if (trafo)
  {
    # transform p-values to Z statistics
    # Find the pvalue columns
    pCols = grep(spaste("^", metaP), colnames(comb));
    if (length(pCols)==0) stop("Internal error: no columns for meta-analysis found. Sorry!");
    setP = comb[, pCols];
    # Caution: I assume here that the returned p-values are two-sided.
    setZ = sign(setStats) * qnorm(setP/2, lower.tail = FALSE);
  } else {
    setZ = setStats;
  }

  colnames(setZ) = spaste("Z.", setNames);
  nObsCols = grep("nPresentSamples", colnames(comb));
  nObs = comb[, nObsCols];

  powers = c(0, 0.5, 1);
  nPowers = 3;

  metaNames = c("equalWeights", "RootDoFWeights", "DoFWeights")
  if (is.null(metaAnalysisWeights)) {
    nMeta = nPowers;
  } else {
    nMeta = nPowers + 1;
    metaNames = c(metaNames, "userWeights");
  }
  metaResults = NULL;
  for (m in 1:nMeta)
  {
    if (m<=nPowers) {
      weights = nObs^powers[m]
    } else
      weights = matrix( metaAnalysisWeights, size$nGenes, nSets, byrow = TRUE);

    metaZ = rowSums( setZ * weights, na.rm = TRUE) / sqrt(rowSums(weights^2, na.rm = TRUE))
    p.meta = 2*pnorm(abs(metaZ), lower.tail = FALSE);
    if (getQvalues)
    {
      q.meta = qvalue.restricted(p.meta);
      meta1 = cbind(metaZ, p.meta, q.meta)
    } else {
      q.meta = NULL;
      meta1 = cbind(metaZ, p.meta);
    }
    colnames(meta1) = spaste(c("Z.", "p.", "q.")[1:ncol(meta1)],
                             metaNames[m]);
    metaResults = cbind(metaResults, meta1);
  }

  # Use rankPvalue to produce yet another meta-analysis

  rankMetaResults = NULL;
  if (useRankPvalue)
  {
    rankPvalueOptions$datS = as.data.frame(setZ);
    if (is.na(match("calculateQvalue", names(rankPvalueOptions))))
      rankPvalueOptions$calculateQvalue = getQvalues;
    for (m in 1:nMeta)
    {
      if (m<=nPowers) {
        weights = nObs^powers[m]
      } else
        weights = matrix( metaAnalysisWeights, size$nGenes, nSets, byrow = TRUE);

      # rankPvalue requires a vector of weights... so compress the weights to a vector.
      # Output a warning if the compression loses information.
      nDifferent = apply(weights, 2, function(x) {length(unique(x)) });
      if (any(nDifferent)>1)
        printFlush(paste("Warning in metaAnalysis: rankPvalue requires compressed weights.\n", 
                         "Some weights may not be entirely accurate."));
      rankPvalueOptions$columnweights = colMeans(weights, na.rm = TRUE);
      rankPvalueOptions$columnweights = rankPvalueOptions$columnweights / sum(rankPvalueOptions$columnweights)
      rp = do.call(rankPvalue, rankPvalueOptions);
      colnames(rp) = spaste(colnames(rp), ".", metaNames[m]);
      rankMetaResults = cbind(rankMetaResults, as.matrix(rp));
    }
  }

  # Put together the output

  out = list(ID = ID,
             metaResults,
             rankMetaResults,
             comb,
             if (trafo) setZ else NULL,
             NULL);   # The last NULL is necessary so the line below works even if nothing else is NULL

  out = as.data.frame(out[ -(which(sapply(out,is.null),arr.ind=TRUE))])

  out;
}


#===============================================================================================
#

