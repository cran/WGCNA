# New multilevel specification of consensus: a hierarchical list.

# The consensus calculation needs to be general enough so it can be used for module merging (consensus
# eigengene network) as well as for calculation of consensus kMEs, and possibly for other stuff as well.

# Consensus specification for a single operation: 
#   - inputs: a list specifying the input of the consensus. This should be general enough to handle
#   blockwise data but also not specific to adjacencies.
#     The inut should be a multiData structure, with each 
#   - calibrationMethod: currently one of "none", "simple quantile", "full quantile"
#   - consensusQuantile
#   - saveCalibratedIndividualTOMs (= FALSE)
#   - calibratedIndividualTOMFilePattern (= "calibratedIndividualTOM-Set%s-Block%b.RData")

# The consensus calculation itself does not need information about correlation type etc; the consensus
# calculation is very general.

# For network analysis applications of the consensus: will also have to keep information about how the
# individual network adjacencies (TOMs) are to be created or were created - correlation type, network type,
# TOM type, correlation options etc.


# So we will keep 2 different types of information around: 
# 1. network construction options. A simple list giving the necessary network construction options.
# 2. ConsensusOptions: This will also be a list giving the options but not holding the actual consensus
# inputs or outputs. 

# outputs of network construction and consensus construction should be held in separate lists.


# Output of individualTOMs: 
#   - adjacency data, a list of blockwiseData instances
#   - block information
#   - network construction options, separately for each adjacency (network construction could be different)

# For consensusTOM:
#   Inputs: 
#      - adjacency data: either from individual TOMs or from other consensus TOMs 
#          . note that constructing a complicated consensus manually (i.e., using consensus TOMs 
#            as inputs to higher-level consensus) is of limited use 
#            since consensus modules would need the full consensus tree anyway.
#      - consensus tree
#   Not really needed: block information
#   outputs:
#      - consensus adjacency data
#      - copy of consensus options
#      - other (diagnostic etc) information

# For consensus modules
#   Inputs: 
#      - undelying expression data
#      - optional consensus TOM
#      - block information
#      - network options for each expression data set
#      - consensus tree


.checkPower = function(power) 
{
    if (any(!is.finite(power)) | (sum(power<1)>0) | (sum(power>50)>0) )  
        stop("power must be between 1 and 50.");
}


#==========================================================================================================
#
# Defining a single consensus operation
#
#==========================================================================================================

newConsensusTree = function(consensusOptions = newConsensusOptions(),
                            inputs,
                            analysisName = NULL)
{
  if (!inherits(consensusOptions, "ConsensusOptions"))
    stop("'consensusOptions' must be of class ConsensusOptions.");
  out = list(consensusOptions = consensusOptions,
       inputs = inputs,
       analysisName = analysisName);
  class(out) = c("ConsensusTree", class(out));
  out;
}

newConsensusOptions = function(
      calibration = c("full quantile", "single quantile", "none"),

      # Simple quantile scaling options
      calibrationQuantile = 0.95,
      sampleForCalibration = TRUE, 
      sampleForCalibrationFactor = 1000,

      # Consensus definition
      consensusQuantile = 0,
      useMean = FALSE,
      setWeights = NULL,
      # Name to prevent clashing of files
      analysisName = "")
{
  calibration = match.arg(calibration);
  if (any(!is.finite(setWeights))) stop("Entries of 'setWeights' must all be finite.");
  if ( (consensusQuantile < 0) | (consensusQuantile > 1) ) 
      stop("'consensusQuantile' must be between 0 and 1.");
  out = list( calibration = calibration,
              calibrationQuantile = calibrationQuantile,
              sampleForCalibration = sampleForCalibration,
              sampleForCalibrationFactor = sampleForCalibrationFactor,
              consensusQuantile = consensusQuantile,
              useMean = useMean,
              setWeights = setWeights,
              analysisName = analysisName);
  class(out) = c("ConsensusOptions", class(out))
  out;
}

newCorrelationOptions = function(
      corType = c("pearson", "bicor"),
      maxPOutliers = 0.05,
      quickCor = 0,
      pearsonFallback = "individual",
      cosineCorrelation = FALSE,
      nThreads = 0,
      corFnc = if (corType=="bicor") "bicor" else "cor",
      corOptions = c(
        list(use = 'p',
             cosine = cosineCorrelation, 
             quick = quickCor,
             nThreads = nThreads), 
        if (corType=="bicor") 
           list(maxPOutliers = maxPOutliers,  
                pearsonFallback = pearsonFallback) else NULL))
{ 
  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) stop("maxPOutliers must be between 0 and 1.");
  if (quickCor < 0) stop("quickCor must be positive.");
  if (nThreads < 0) stop("nThreads must be positive.");
  corType = match.arg(corType);
  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) stop("maxPOutliers must be between 0 and 1.");
  if (quickCor < 0) stop("quickCor must be positive.");
  fallback = pmatch(pearsonFallback, .pearsonFallbacks)
  if (is.na(fallback))
      stop(paste("Unrecognized 'pearsonFallback'. Recognized values are (unique abbreviations of)\n",
           paste(.pearsonFallbacks, collapse = ", ")))
  out = list(
    corType = corType,
    maxPOutliers = maxPOutliers,
    quickCor = quickCor,
    pearsonFallback = pearsonFallback,
    pearsonFallback.code = match(pearsonFallback, .pearsonFallbacks),
    cosineCorrelation = cosineCorrelation,
    corFnc = corFnc,
    corOptions = corOptions,
    corType.code = match(corType, .corTypes),
    nThreads = nThreads);
  class(out) = c("CorrelationOptions", class(out));
  out;
}

newNetworkOptions = function(
    correlationOptions = newCorrelationOptions(),
 
    # Adjacency options
    replaceMissingAdjacencies = TRUE,
    power = 6,
    networkType = c("signed hybrid", "signed", "unsigned"),
    checkPower = TRUE,

    # Topological overlap options
    TOMType = c("signed", "unsigned", "none"),
    TOMDenom = c("mean", "min"),
    suppressTOMForZeroAdjacencies = FALSE,

    # Internal behavior options
    useInternalMatrixAlgebra = FALSE)
{
  if (checkPower) .checkPower(power);
  networkType = match.arg(networkType);
  TOMType = match.arg(TOMType);
  TOMDenom = match.arg(TOMDenom);
  out = c(correlationOptions,
      list(replaceMissingAdjacencies = replaceMissingAdjacencies,
           power = power,
           networkType = networkType,
           TOMType = TOMType,
           TOMDenom = TOMDenom,
           networkType.code = match(networkType, .networkTypes),
           TOMType.code = match(TOMType, .TOMTypes),
           TOMDenom.code = match(TOMDenom, .TOMDenoms),
           suppressTOMForZeroAdjacencies = suppressTOMForZeroAdjacencies,
           useInternalMatrixAlgebra = useInternalMatrixAlgebra));
  class(out) = c("NetworkOptions", class(correlationOptions));
  out;
}

#====================================================================================================
#
# cor, network, and consensus calculations
#
#====================================================================================================

.corCalculation = function(x, y = NULL, weights.x = NULL, weights.y = NULL, correlationOptions)
{
  if (!inherits(correlationOptions, "CorrelationOptions"))
    stop(".corCalculation: 'correlationOptions' does not have correct type.");
  do.call(match.fun(correlationOptions$corFnc), 
          c(list(x = x, y= y, weights.x = weights.x, weights.y = weights.y), 
            correlationOptions$corOptions));
}


# network calculation. Returns the resulting topological overlap or adjacency matrix.

.networkCalculation = function(data, networkOptions, weights = NULL, 
                verbose = 1, indent = 0)
{
  if (is.data.frame(data)) data = as.matrix(data);
  if (!inherits(networkOptions, "NetworkOptions"))
    stop(".networkCalculation: 'networkOptions' does not have correct type.");
 
   .checkAndScaleWeights(weights, data, scaleByMax = FALSE);

   callVerb = max(0, verbose - 1); callInd = indent + 2;
   CcorType = networkOptions$corType.code - 1;
   CnetworkType = networkOptions$networkType.code - 1;
   CTOMType = networkOptions$TOMType.code -1;
   CTOMDenom = networkOptions$TOMDenom.code -1;

   warn = as.integer(0);
   # For now return the full matrix; eventually we may return just the dissimilarity.
   # To make full use of the lower triagle space saving we'd have to modify the underlying C code
   # dramatically, otherwise it will still need to allocate the full matrix for the matrix multiplication.
   .Call("tomSimilarity_call", data, weights,
          as.integer(CcorType), as.integer(CnetworkType), 
          as.double(networkOptions$power), as.integer(CTOMType),
          as.integer(CTOMDenom),
          as.double(networkOptions$maxPOutliers),
          as.double(networkOptions$quickCor),
          as.integer(networkOptions$pearsonFallback.code),
          as.integer(networkOptions$cosineCorrelation),
          as.integer(networkOptions$replaceMissingAdjacencies),
          as.integer(networkOptions$suppressTOMForZeroAdjacencies),
          as.integer(networkOptions$useInternalMatrixAlgebra),
          warn, as.integer(min(1, networkOptions$nThreads)),
          as.integer(callVerb), as.integer(callInd), PACKAGE = "WGCNA");
}

# the following is contained in consensusOptions:
#  out = list( calibration = calibration,
#              calibrationQuantile = calibrationQuantile,
#              sampleForCalibration = sampleForCalibration,
#              sampleForCalibrationFactor = sampleForCalibrationFactor,
#              consensusQuantile = consensusQuantile,
#              useMean = useMean,
#              setWeights = setWeights);


#==============================================================================================
#
# general utility functions
#
#==============================================================================================

# Try to guess whether disk cache should be used
# this should work for both multiData as well as for simple lists of arrays.

.useDiskCache = function(multiExpr, blocks = NULL, chunkSize = NULL,
                         nSets = if (isMultiData(multiExpr)) length(multiExpr) else 1,
                         nGenes = checkSets(multiExpr)$nGenes)
{
  if (is.null(chunkSize)) chunkSize = as.integer(.largestBlockSize/(2*nSets))

  if (length(blocks) == 0) {
    blockLengths = nGenes;
  } else 
    blockLengths = as.numeric(table(blocks));

  max(blockLengths) > chunkSize;
}

.dimensions = function(x)
{
   if (is.null(dim(x))) return(length(x))
   return(dim(x))
}

.shiftList = function(c0, lst)
{
  if (length(lst)==0) return(list(c0));
  ll = length(lst)
  out = list();
  out[[1]] = c0;
  for (i in 1:ll)
    out[[i+1]] = lst[[i]];
  out;
}

.checkListDimConsistencyAndGetDimnames = function(dataList)
{
  nPars = length(dataList)
  dn = NULL
  for (p in 1:nPars)
  {
      if (mode(dataList[[p]])!="numeric")
          stop(paste("Argument number", p, " is not numeric."))
       if (p==1) {
          dim = .dimensions(dataList[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(dataList[[p]]), dim)))
             stop("Argument dimensions are not consistent.");
       }
       if (prod(dim)==0) stop(paste("Argument has zero dimension."));
       if (is.null(dn)) dn = dimnames(dataList[[p]]);
   }
   dn;
}

.mtd.checkDimConsistencyAndGetDimnames = function(mtd)
{
  nPars = length(mtd)
  dn = NULL
  for (p in 1:nPars)
  {
      if (mode(mtd[[p]]$data)!="numeric")
          stop(paste("Argument number", p, " is not numeric."))
       if (p==1) {
          dim = .dimensions(mtd[[p]]$data)
       } else {
          if (!isTRUE(all.equal(.dimensions(mtd[[p]]$data), dim)))
             stop("Argument dimensions are not consistent.");
       }
       if (prod(dim)==0) stop(paste("Argument has zero dimension."));
       if (is.null(dn)) dn = dimnames(mtd[[p]]$data);
   }
   dn;
}



#==============================================================================================
#
# general utility functions for working with disk cache
#
#==============================================================================================

.saveChunks = function(data, chunkSize, fileBase, cacheDir, fileNames = NULL)
{
   ld = length(data);
   nChunks = ceiling(ld/chunkSize);
   if (is.null(fileNames))
   {
     if (length(fileBase)!=1)
       stop("Internal error: length of 'fileBase' must be 1.");

     fileNames = rep("", nChunks);
     x = 1;
     for (c in 1:nChunks)
     {
       fileNames[c] = tempfile(pattern = fileBase, tmpdir = cacheDir);
       # This effectively reserves the file name
       save(x, file = fileNames[c]);
     }
   } else {
     if (length(fileNames)!=nChunks)
       stop("Internal error: length of 'fileNames' must equal the number of chunks.");
   }

   chunkLengths = rep(0, nChunks);
   start = 1;
   for (c in 1:nChunks)
   {
     end = min(start + chunkSize-1, ld);
     chunkLengths[c] = end - start + 1;
     temp = data[start:end];
     save(temp, file = fileNames[c]);
     start = end + 1;
   }
   rm(temp); #gc();
   list(files = fileNames, chunkLengths = chunkLengths);
}

.loadAsList = function(file)
{
  env = new.env();
  load(file = file, envir = env);
  as.list(env);
};


.loadObject = function(file, name = NULL, size = NULL)
{
  x = .loadAsList(file);
  if (!is.null(name) && (names(x)[1]!=name))
    stop("File ", file, " does not contain object '", name, "'.")

  obj = x[[1]];
  if (!is.null(size) && (length(obj)!=size))
    stop("Object '", name, "' does not have the correct length.");
  obj;
}

.vector2dist = function(x)
{
  n = length(x);
  n1 = (1 + sqrt(1 + 8*n))/2
  if (floor(n1)!=n1) stop("Input length not consistent with a distance structure.");
  attributes(x) = list(Size = as.integer(n1), Diag = FALSE, Upper = FALSE);
  class(x) = "dist";
  x;
}

.emptyDist = function(nObjects, fill = 0)
{
  n = (nObjects * (nObjects-1))/2;
  .vector2dist(rep(fill, n));
}

.checkAndDelete = function(files)
{
  if (length(files)>0) lapply(as.character(files), function(file) if (file.exists(file)) file.remove(file));
  NULL;
}

.qorder = function(data)
{
  data = as.numeric(data);
  .Call("qorder", data, PACKAGE = "WGCNA")
}

# Actual consensus calculation distilled into one function. data is assumed to have sets in columns
# and samples/observations/whatever in rows. setWeightMat should be a matrix of dimensions (nSets, 1) 
# and be normalized to sum=1.

.consensusCalculation.base = function(data, useMean, setWeightMat, consensusQuantile)
{
  nSets = ncol(data);
  if (nSets==1)
  {
    out.list = list(consensus = data);
    if (consensusQuantile==0) out.list$originCount = c(`1`=nrow(data))
  } else {
    if (useMean)
    {
       if (any(is.na(data)))
       {
         finiteMat = 1-is.na(data);
         data[is.na(data)] = 0;
         out = data %*% setWeightMat / finiteMat%*%setWeightMat;
       } else {
         out = data %*% setWeightMat;
       }
       out.list = list(consensus = out);
    } else if (consensusQuantile == 0) 
    {
        #min = rep(0, nrow(data));
        #which = rep(0, nrow(data));
        #whichmin = .C("minWhichMin_row", as.double(data),
        #              as.integer(nrow(data)), as.integer(ncol(data)),
        #              as.double(min), as.double(which), PACKAGE = "WGCNA");
        #min = whichmin[[4]];
        #which = whichmin[[5]] + 1;
        #rm(whichmin);
        whichmin = .Call("minWhich_call", data, 1L, PACKAGE = "WGCNA");
        out.list = list(consensus = whichmin$min, originCount = table(as.integer(whichmin$which)));
    } else {
       out.list = list(consensus = rowQuantileC(data, p = consensusQuantile));
    }
  }
  out.list;
}

.consensusCalculation.base.FromList = function(dataList, useMean, setWeights, consensusQuantile)
{
  nSets = length(dataList);
  if (nSets==1)
  {
    out.list = list(consensus = dataList[[1]]);
    if (consensusQuantile == 0) out.list$originCount = c(`1`=length(dataList[[1]]))
  } else {
    if (useMean)
    {
       out.list = list(consensus = pmean(dataList, setWeights));
    } else if (consensusQuantile == 0)
    {
        whichmin = pminWhich.fromList(dataList);
        min = whichmin$min;
        which = whichmin$which;
        out.list = list(consensus = min, originCount = table(as.integer(which)));
    } else {
       out.list = list(consensus = pquantile.fromList(dataList, prob = consensusQuantile));
    }
  }
  out.list;
}

#==============================================================================================
#
# utility functions for working with multiple blockwise adjacencies.
#
#==============================================================================================

newBlockInformation = function(
    blocks,
    goodSamplesAndGenes)
{
  blockGenes = tapply(1:length(blocks), blocks, identity)
  names(blockGenes) = sort(unique(blocks));
  nGGenes = sum(goodSamplesAndGenes$goodGenes);
  gBlocks = blocks[goodSamplesAndGenes$goodGenes];
  out = list(blocks = blocks,
             blockGenes = blockGenes, 
             goodSamplesAndGenes = goodSamplesAndGenes,
             nGGenes = nGGenes,
             gBlocks = gBlocks);
  class(out) = c("BlockInformation", class(out));
  out;
}
  
#=======================================================================================================
#
# individualTOMs
# This is essentially a re-badged blockwiseIndividualTOMs with a different input and ouptut format.
#
#=======================================================================================================

# The following is contained in networkOptions:

    #corType = corType,
    #maxPOutliers = maxPOutliers,
    #quickCor = quickCor,
    #pearsonFallback = pearsonFallback,
    #cosineCorrelation = cosineCorrelation,
    #corFnc = corFnc,
    #corOptions = corOptions,
    #corType.code = match(corType, .corTypes),

    # Adjacency options
    #replaceMissingAdjacencies = TRUE,
    #power = 6,
    #networkType = c("signed hybrid", "signed", "unsigned"),

    # Topological overlap options

    #TOMType = c("signed", "unsigned"),
    #TOMDenom = c("mean", "min"))


individualTOMs = function(
   multiExpr,
   multiWeights = NULL,

   multiExpr.imputed = NULL,  ## Optional, useful for pre-clustering if preclustering is needed.

   # Data checking options

   checkMissingData = TRUE,

   # Blocking options

   blocks = NULL, 
   maxBlockSize = 5000, 
   blockSizePenaltyPower = 5,
   nPreclusteringCenters = NULL,
   randomSeed = 12345,

   # Network construction options. This can be a single object of class NetworkOptions, or a multiData
   # structure of NetworkOptions objects, one per element of multiExpr.

   networkOptions,

   # Save individual TOMs? This is equivalent to using external = TRUE in blockwiseData

   saveTOMs = TRUE,
   individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",

   # Behaviour options
   collectGarbage = TRUE,
   verbose = 2, indent = 0)
{
  spaces = indentSpaces(indent);

  dataSize = checkSets(multiExpr, checkStructure = TRUE);
  if (dataSize$structureOK)
  {
    nSets = dataSize$nSets;
    nGenes = dataSize$nGenes;
    multiFormat = TRUE;
  } else {
    multiExpr = multiData(multiExpr);
    if (!is.null(multiWeights)) multiWeights = multiData(multiWeights);
    nSets = dataSize$nSets;
    nGenes = dataSize$nGenes;
    multiFormat = FALSE;
  }

  .checkAndScaleMultiWeights(multiWeights, multiExpr, scaleByMax = FALSE);

  if (inherits(networkOptions, "NetworkOptions"))
  {
    networkOptions = list2multiData(.listRep(networkOptions, nSets));
  } else {
    if (!all(mtd.apply(networkOptions, inherits, "NetworkOptions", mdaSimplify = TRUE)))
       stop("'networkOptions' must be of class 'NetworkOptions' or a multiData structure\n",
            "   of objects of the class.\n",
            "   See newNetworkOptions for creating valid network options.");
  }

  if (is.null(names(multiExpr)))
    names(multiExpr) = spaste("Set", 1:nSets);

  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed);
    } 
    set.seed(randomSeed);
  }

  #if (maxBlockSize >= floor(sqrt(2^31)) )
  #  stop("'maxBlockSize must be less than ", floor(sqrt(2^31)), ". Please decrease it and try again.")

  if (!is.null(blocks) && (length(blocks)!=nGenes))
    stop("Input error: length of 'blocks' must equal number of genes in 'multiExpr'.");

  if (verbose>0) 
     printFlush(paste(spaces, "Calculating topological overlaps block-wise from all genes"));

  nSamples = dataSize$nSamples;
  
  # Check data for genes and samples that have too many missing values

  # Check that multiExpr has valid (mtd)column names. If column names are missing, generate them.

  colIDs = mtd.colnames(multiExpr);
  if (is.null(colIDs)) colIDs = c(1:dataSize$nGenes);

  if (checkMissingData)
  {
    gsg = goodSamplesGenesMS(multiExpr, multiWeights = multiWeights, 
                             verbose = verbose - 1, indent = indent + 1)
    if (!gsg$allOK)
    {
      multiExpr = mtd.subset(multiExpr, gsg$goodSamples, gsg$goodGenes);
      if (!is.null(multiWeights)) 
        multiWeights = mtd.subset(multiWeights, gsg$goodSamples, gsg$goodGenes);
      if (!is.null(multiExpr.imputed)) 
        multiExpr.imputed = mtd.subset(multiExpr.imputed, gsg$goodSamples, gsg$goodGenes);
      colIDs = colIDs[gsg$goodGenes];
    }
  } else {
    gsg = list(goodGenes = rep(TRUE, nGenes), goodSamples = lapply(nSamples, function(n) rep(TRUE, n)));
    gsg$allOK = TRUE;
  }
  nGGenes = sum(gsg$goodGenes);
  nGSamples = rep(0, nSets);
  for (set in 1:nSets) nGSamples[set] = sum(gsg$goodSamples[[set]]);

  if (is.null(blocks))
  {
    if (nGGenes > maxBlockSize)
    {
      if (verbose>1) printFlush(paste(spaces, "....pre-clustering genes to determine blocks.."));
      clustering = consensusProjectiveKMeans(
                       if (!is.null(multiExpr.imputed)) multiExpr.imputed else multiExpr, 
                       preferredSize = maxBlockSize,
                       sizePenaltyPower = blockSizePenaltyPower, checkData = FALSE,
                       nCenters = nPreclusteringCenters,
                       verbose = verbose-2, indent = indent + 1);
      gBlocks = .orderLabelsBySize(clustering$clusters);
    } else 
      gBlocks = rep(1, nGGenes);
    blocks = rep(NA, nGenes);
    blocks[gsg$goodGenes] = gBlocks;
  } else {
    gBlocks = blocks[gsg$goodGenes];
  }

  blockLevels = as.numeric(levels(factor(gBlocks)));
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels);

  if (any(blockSizes > sqrt(2^31)-1))
    printFlush(spaste(spaces,
            "Found block(s) with size(s) larger than limit of 'int' indexing.\n",
            spaces, " Support for such large blocks is experimental; please report\n",
            spaces, " any issues to Peter.Langfelder@gmail.com."));

  # check file names for uniqueness

  actualFileNames = matrix("", nSets, nBlocks);
  if (saveTOMs) 
  {
    for (set in 1:nSets) for (b in 1:nBlocks)
      actualFileNames[set, b] = .processFileName(individualTOMFileNames, set, names(multiExpr), b);
  
    rownames(actualFileNames) = spaste("Set.", c(1:nSets));
    colnames(actualFileNames) = spaste("Block.", c(1:nBlocks));
    if (length(unique(as.vector(actualFileNames))) < nSets * nBlocks) 
    {
      printFlush("Error: File names for (some) set/block combinations are not unique:");
      print(actualFileNames);
      stop("File names must be unique.");
    }
  }
  
  # Initialize various variables

  blockGenes = list();
  blockNo = 1;
  setTomDS = replicate(nSets, list());
  # Here's where the analysis starts
  for (blockNo in 1:nBlocks)
  {
    if (verbose>1 && nBlocks > 1) printFlush(paste(spaces, "..Working on block", blockNo, "."));
    # Select the block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]];
    #nBlockGenes = length(block);
    #blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks==blockLevels[blockNo]];
    #errorOccurred = FALSE;

    # For each set: calculate and save TOM

    for (set in 1:nSets)
    {
      if (verbose>2) printFlush(paste(spaces, "....Working on set", set))
      selExpr = as.matrix(multiExpr[[set]]$data[, block]);
      if (!is.null(multiWeights)) {
         selWeights = as.matrix(multiWeights[[set]]$data[, block]) 
      } else
         selWeights = NULL;

      tomDS = as.dist(.networkCalculation(selExpr, networkOptions[[set]]$data, weights = selWeights, 
                      verbose = verbose -2, indent = indent+2));

      setTomDS[[set]]$data = addBlockToBlockwiseData(if (blockNo==1) NULL else setTomDS[[set]]$data, 
                            external = saveTOMs,
                            blockData = tomDS, 
                            blockFile = actualFileNames[set, blockNo],
                            recordAttributes = TRUE,
                            metaData = list(IDs = colIDs[block]))
    }
    if (collectGarbage) { rm(tomDS); gc(); }
  }

  names(setTomDS) = names(multiExpr);

  if (!multiFormat)
  {
    gsg$goodSamples = gsg$goodSamples[[1]];
    setTomDS = setTomDS[[1]]$data;
  }

  blockInformation = newBlockInformation(blocks, gsg);

  list(blockwiseAdjacencies = setTomDS,
       setNames = names(multiExpr),
       nSets = length(multiExpr),
       blockInfo = blockInformation,
       networkOptions = networkOptions
       )
}


#=====================================================================================================
#
# hierarchical consensus TOM
#
#=====================================================================================================

hierarchicalConsensusTOM = function(
      # Supply either ...
      # ... information needed to calculate individual TOMs

      multiExpr,
      multiWeights = NULL,

      # Data checking options
      checkMissingData = TRUE,

      # Blocking options
      blocks = NULL,
      maxBlockSize = 20000,
      blockSizePenaltyPower = 5,
      nPreclusteringCenters = NULL,
      randomSeed = 12345,

      # Network construction options. This can be a single object of class NetworkOptions, or a multiData
      # structure of NetworkOptions objects, one per element of multiExpr.

      networkOptions,

      # Save individual TOMs?

      keepIndividualTOMs = TRUE,
      individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",

      # ... or information about individual (more precisely, input) TOMs

      individualTOMInfo = NULL,

      # Consensus calculation options 
      consensusTree,

      useBlocks = NULL,

      # Save calibrated TOMs?
      saveCalibratedIndividualTOMs = FALSE,
      calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set%s-Block%b.RData",

      # Return options
      saveConsensusTOM = TRUE,
      consensusTOMFilePattern = "consensusTOM-%a-Block%b.RData",
      getCalibrationSamples = FALSE,

      # Return the intermediate results as well?  
      keepIntermediateResults = saveConsensusTOM,

      # Internal handling of TOMs
      useDiskCache = NULL, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      # Behavior
      collectGarbage = TRUE,
      verbose = 1,
      indent = 0)
{
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed);
    } 
    set.seed(randomSeed);
  }

  localIndividualTOMCalculation = is.null(individualTOMInfo);


  if (is.null(individualTOMInfo))
  {
    if (missing(multiExpr)) stop("Either 'individualTOMInfo' or 'multiExpr' must be given.");
    if (is.null(useDiskCache)) useDiskCache = .useDiskCache(multiExpr, blocks, chunkSize);
    time = system.time({individualTOMInfo = individualTOMs(multiExpr = multiExpr, 
                         multiWeights = multiWeights,
                         checkMissingData = checkMissingData,

                         blocks = blocks,
                         maxBlockSize = maxBlockSize,
                         blockSizePenaltyPower = blockSizePenaltyPower,
                         nPreclusteringCenters = nPreclusteringCenters,
                         randomSeed = NULL,

                         networkOptions = networkOptions,

                         saveTOMs = useDiskCache | keepIndividualTOMs,
                         individualTOMFileNames = individualTOMFileNames,

                         collectGarbage = collectGarbage,
                         verbose = verbose, indent = indent);});
    if (verbose > 1) { printFlush("Timimg for individual TOMs:"); print(time); }
  }
  consensus = hierarchicalConsensusCalculation(individualTOMInfo$blockwiseAdjacencies,
           consensusTree,
           level = 1,
           useBlocks = useBlocks,
           randomSeed = NULL,
           saveCalibratedIndividualData = saveCalibratedIndividualTOMs,
           calibratedIndividualDataFilePattern = calibratedIndividualTOMFilePattern,
           saveConsensusData = saveConsensusTOM,
           consensusDataFileNames = consensusTOMFilePattern,
           getCalibrationSamples= getCalibrationSamples,
           keepIntermediateResults = keepIntermediateResults,
           useDiskCache = useDiskCache, 
           chunkSize = chunkSize,
           cacheDir = cacheDir,
           cacheBase = cacheBase,
           collectGarbage = collectGarbage,
           verbose = verbose,
           indent = indent);
           
  if (localIndividualTOMCalculation)
  {
     if (!keepIndividualTOMs)
     {
        # individual TOMs are contained in the individualTOMInfo list; remove them.
        mtd.apply(individualTOMInfo$blockwiseAdjacencies, BD.checkAndDeleteFiles);
        individualTOMInfo$blockwiseAdjacencies = NULL
     }
  } 
        
  c( consensus,
     list(individualTOMInfo = individualTOMInfo,
          consensusTree = consensusTree)
   );
}


#========================================================================================================
#
# Merge consensusTOMInfo lists
#
#========================================================================================================
#
# Caution: at present the function does not check that the inputs are compatible and compatible with the
# supplied blocks.

.mergeConsensusTOMInformationLists = function(blockInformation, consensusTOMInfoList)
{

  blocks = blockInformation$blocks;

  blockLevels = unique(blocks);
  nBlocks = length(blockLevels);

  out = consensusTOMInfoList[[1]];

  # Merging consensus information
  out$consensusData = do.call(mergeBlockwiseData, lapply(consensusTOMInfoList, getElement, "consensusData"));
  if (out$saveCalibratedIndividualData)
  {
    out$calibratedIndividualData = lapply(1:out$nSets, function(set)
           do.call(mergeBlockwiseData, 
                lapply(consensusTOMInfoList, function(ct) ct$calibratedIndividualData[[set]])));
  }

  if (!is.null(out$calibrationSamples))
  {
    out$calibrationSamples = do.call(c, lapply(consensusTOMInfoList, getElement, "calibrationSamples"));
  }
  out$originCount = rowSums(do.call(cbind, lapply(consensusTOMInfoList, getElement, "originCount")))

  # Merging information in individualTOMInfo

  out$individualTOMInfo$blockwiseAdjacencies = lapply(1:out$nSets, function(set)
      list(data = do.call(mergeBlockwiseData, lapply(consensusTOMInfoList, function(ct)
                     ct$individualTOMInfo$blockwiseAdjacencies[[set]]$data)))); 

  out$individualTOMInfo$blockInfo = blockInformation;

  out;
}

#========================================================================================================
#
# consensusTOM: old, single-layer consensus.
#
#========================================================================================================

consensusTOM = function(
      # Supply either ...
      # ... information needed to calculate individual TOMs

      multiExpr,

      # Data checking options
      checkMissingData = TRUE,

      # Blocking options
      blocks = NULL,
      maxBlockSize = 5000,
      blockSizePenaltyPower = 5,
      nPreclusteringCenters = NULL,
      randomSeed = 12345,

      # Network construction arguments: correlation options

      corType = "pearson",
      maxPOutliers = 1,
      quickCor = 0,
      pearsonFallback = "individual",
      cosineCorrelation = FALSE,
      replaceMissingAdjacencies = FALSE,

      # Adjacency function options

      power = 6,
      networkType = "unsigned",
      checkPower = TRUE,

      # Topological overlap options

      TOMType = "unsigned",
      TOMDenom = "min",

      # Save individual TOMs?

      saveIndividualTOMs = TRUE,
      individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",

      # ... or individual TOM information

      individualTOMInfo = NULL,
      useIndivTOMSubset = NULL,

   ##### Consensus calculation options 

      useBlocks = NULL,

      networkCalibration = c("single quantile", "full quantile", "none"),

      # Save calibrated TOMs?
      saveCalibratedIndividualTOMs = FALSE,
      calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set%s-Block%b.RData",

      # Simple quantile scaling options
      calibrationQuantile = 0.95,
      sampleForCalibration = TRUE, sampleForCalibrationFactor = 1000,
      getNetworkCalibrationSamples = FALSE,

      # Consensus definition
      consensusQuantile = 0,
      useMean = FALSE,
      setWeights = NULL,

      # Return options
      saveConsensusTOMs = TRUE,
      consensusTOMFilePattern = "consensusTOM-Block%b.RData",
      returnTOMs = FALSE,

      # Internal handling of TOMs
      useDiskCache = NULL, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      nThreads = 1,

      # Diagnostic messages
      verbose = 1,
      indent = 0)
{
  spaces = indentSpaces(indent);
  networkCalibration = match.arg(networkCalibration);

  seedSaved = FALSE;
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed);
    } 
    set.seed(randomSeed);
  }

  if (any(!is.finite(setWeights))) stop("Entries of 'setWeights' must all be finite.");

  localIndividualTOMCalculation = is.null(individualTOMInfo);
  if (is.null(individualTOMInfo))
  {
    if (missing(multiExpr)) stop("Either 'individualTOMInfo' or 'multiExpr' must be given.");

    dataSize = checkSets(multiExpr);
    nSets.all = dataSize$nSets;
    nGenes = dataSize$nGenes;

    if (is.null(useDiskCache)) useDiskCache = .useDiskCache(multiExpr, blocks, chunkSize);
    if (length(power)!=1)
    {
      if (length(power)!=nSets.all)
        stop("Invalid arguments: Length of 'power' must equal number of sets given in 'multiExpr'.");
    } else {
      power = rep(power, nSets.all);
    }

    if ( (consensusQuantile < 0) | (consensusQuantile > 1) ) 
      stop("'consensusQuantile' must be between 0 and 1.");

    time = system.time({individualTOMInfo = blockwiseIndividualTOMs(multiExpr = multiExpr, 
                         checkMissingData = checkMissingData,
                         blocks = blocks,
                         maxBlockSize = maxBlockSize,
                         blockSizePenaltyPower = blockSizePenaltyPower,
                         nPreclusteringCenters = nPreclusteringCenters,
                         randomSeed = NULL,
                         corType = corType,
                         maxPOutliers = maxPOutliers,
                         quickCor = quickCor,
                         pearsonFallback = pearsonFallback,
                         cosineCorrelation = cosineCorrelation,
                         replaceMissingAdjacencies = replaceMissingAdjacencies,
                         power = power,
                         networkType = networkType, 
                         TOMType = TOMType,
                         TOMDenom = TOMDenom,
                         saveTOMs = useDiskCache | saveIndividualTOMs,
                         individualTOMFileNames = individualTOMFileNames,
                         nThreads = nThreads,
                         verbose = verbose, indent = indent);});
    if (verbose > 1) { printFlush("Timimg for individual TOMs:"); print(time); }

    if (!saveIndividualTOMs & useDiskCache) 
       on.exit(.checkAndDelete(individualTOMInfo$actualTOMFileNames), add = TRUE)

  } else {
    nSets.all = if (individualTOMInfo$saveTOMs)
             nrow(individualTOMInfo$actualTOMFileNames) else ncol(individualTOMInfo$TOMSimilarities[[1]]);
    nGenes = length(individualTOMInfo$blocks);
    if (is.null(useDiskCache)) useDiskCache = .useDiskCache(NULL, blocks, chunkSize, 
               nSets = nSets.all, nGenes = nGenes);

  }
  
  nGoodGenes = length(individualTOMInfo$gBlocks);

  if (is.null(setWeights)) setWeights = rep(1, nSets.all);
  if (length(setWeights)!=nSets.all)
    stop("Length of 'setWeights' must equal the number of sets.");

  setWeightMat = as.matrix(setWeights/sum(setWeights));

  if (is.null(useIndivTOMSubset))
  {  
    if (individualTOMInfo$nSets != nSets.all)
      stop(paste("Number of sets in individualTOMInfo and in multiExpr do not agree.\n",
                 "  To use a subset of individualTOMInfo, set useIndivTOMSubset appropriately."));

    useIndivTOMSubset = c(1:nSets.all);
  }

  nSets = length(useIndivTOMSubset);


  if (length(unique(useIndivTOMSubset))!=nSets)
    stop("Entries of 'useIndivTOMSubset' must be unique");

  if (any(useIndivTOMSubset<1) | any(useIndivTOMSubset>individualTOMInfo$nSets))
    stop("All entries of 'useIndivTOMSubset' must be between 1 and the number of sets in individualTOMInfo");

  # if ( (minKMEtoJoin >1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.");

  gsg = individualTOMInfo$goodSamplesAndGenes;

  # Restrict gsg to used sets

  gsg$goodSamples = gsg$goodSamples[useIndivTOMSubset];

  if (is.null(chunkSize)) chunkSize = as.integer(.largestBlockSize/(2*nSets))

  # Initialize various variables

  if (getNetworkCalibrationSamples)
  {
    if (!sampleForCalibration)
      stop(paste("Incompatible input options: networkCalibrationSamples can only be returned", 
                 "if sampleForCalibration is TRUE."));
    networkCalibrationSamples = list();
  }

  blockLevels = sort(unique(individualTOMInfo$gBlocks));
  nBlocks = length(blockLevels);

  if (is.null(useBlocks)) useBlocks = blockLevels;

  useBlockIndex = match(useBlocks, blockLevels);

  if (!all(useBlocks %in% blockLevels))
    stop("All entries of 'useBlocks' must be valid block levels.");

  if (any(duplicated(useBlocks)))
    stop("Entries of 'useBlocks' must be unique.");

  nUseBlocks = length(useBlocks);
  if (nUseBlocks==0)
    stop("'useBlocks' cannot be non-NULL and empty at the same time.");

  consensusTOM.out = list();

  TOMFiles = rep("", nUseBlocks);
  originCount = rep(0, nSets);

  calibratedIndividualTOMFileNames = NULL;
  if (saveCalibratedIndividualTOMs)
  {
    calibratedIndividualTOMFileNames = matrix("", nSets, nBlocks);
    for (set in 1:nSets) for (b in 1:nBlocks)
      calibratedIndividualTOMFileNames[set, b] = .processFileName(calibratedIndividualTOMFilePattern, 
                                 setNumber = set, blockNumber = b);

  }
  gc();

  # Here's where the analysis starts

  for (blockIndex in 1:nUseBlocks)
  {
    blockNo = useBlockIndex[blockIndex];

    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."));
    # Select block genes
    block = c(1:nGoodGenes)[individualTOMInfo$gBlocks==blockLevels[blockNo]];
    nBlockGenes = length(block);
    # blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks==blockLevels[blockNo]];
    scaleQuant = rep(1, nSets);
    scalePowers = rep(1, nSets);

    # Set up file names or memory space to hold the set TOMs
    if (useDiskCache)
    {
      nChunks = ceiling(nBlockGenes * (nBlockGenes-1)/2/chunkSize);
      chunkFileNames = array("", dim = c(nChunks, nSets));
      on.exit(.checkAndDelete(chunkFileNames), add = TRUE);
    } else nChunks = 1;

    if (nChunks==1) useDiskCache = FALSE;
    if (!useDiskCache)
    {
      # Note: setTomDS will contained the scaled set TOM matrices.
      setTomDS = array(0, dim = c(nBlockGenes*(nBlockGenes-1)/2, nSets));
    } 

    # create an empty consTomDS distance structure.

    consTomDS = .emptyDist(nBlockGenes);

    # sample entry indices from the distance structure for TOM scaling, if requested

    if (networkCalibration=="single quantile" && sampleForCalibration)
    {
      qx = min(calibrationQuantile, 1-calibrationQuantile);
      nScGenes = min(sampleForCalibrationFactor * 1/qx, length(consTomDS));
      nTOMEntries = length(consTomDS)
      scaleSample = sample(nTOMEntries, nScGenes);
      if (getNetworkCalibrationSamples)
        networkCalibrationSamples[[blockIndex]] = list(sampleIndex = scaleSample,
                                            TOMSamples = matrix(NA, nScGenes, nSets));
    }
    if (networkCalibration %in% c("single quantile", "none"))
    {
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(paste(spaces, "....Working on set", useIndivTOMSubset[set]))
        if (individualTOMInfo$saveTOMs)
        {
           tomDS = .loadObject(individualTOMInfo$ actualTOMFileNames[useIndivTOMSubset[set], blockNo],
                               name = "tomDS", size = nBlockGenes*(nBlockGenes-1)/2);
        } else {
          tomDS = consTomDS;
          tomDS[] = individualTOMInfo$TOMSimilarities[[blockNo]] [, useIndivTOMSubset[set]]
        }
        
        if (networkCalibration=="single quantile")
        {
          # Scale TOMs so that calibrationQuantile agree in each set
          if (sampleForCalibration)
          {
            if (getNetworkCalibrationSamples)
            { 
              networkCalibrationSamples[[blockIndex]]$TOMSamples[, set] = tomDS[scaleSample];
              scaleQuant[set] = quantile(networkCalibrationSamples[[blockIndex]]$TOMSamples[, set], 
                                         probs = calibrationQuantile, type = 8);
            } else {
              scaleQuant[set] = quantile(tomDS[scaleSample], probs = calibrationQuantile, type = 8);
            }
          } else
            scaleQuant[set] = quantile(x = tomDS, probs = calibrationQuantile, type = 8);
          if (set>1) 
          {
             scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
             tomDS = tomDS^scalePowers[set];
          }
          if (saveCalibratedIndividualTOMs)
             save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo]);
        } 

        # Save the calculated TOM either to disk in chunks or to memory.
      
        if (useDiskCache)
        {
          if (verbose > 3) printFlush(paste(spaces, "......saving TOM similarity to disk cache.."));
          sc = .saveChunks(tomDS, chunkSize, cacheBase, cacheDir = cacheDir);
          chunkFileNames[, set] = sc$files;
          chunkLengths = sc$chunkLengths;
        } else {
          setTomDS[, set] = tomDS[];
        }
        rm(tomDS); gc();
      }
    } else if (networkCalibration=="full quantile")
    {
      # Step 1: load each TOM, get order, split TOM into chunks according to order, and save.
      if (verbose>1) printFlush(spaste(spaces, "..working on quantile normalization"))
      if (useDiskCache)
      {
        orderFiles = rep("", nSets);
        on.exit(.checkAndDelete(orderFiles),add = TRUE);
      }
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(paste(spaces, "....Working on set", useIndivTOMSubset[set]))
        if (individualTOMInfo$saveTOMs)
        {
           tomDS = .loadObject(individualTOMInfo$ actualTOMFileNames[useIndivTOMSubset[set], blockNo],
                               name = "tomDS", size = nBlockGenes*(nBlockGenes-1)/2);
        } else {
          tomDS = consTomDS;
          tomDS[] = individualTOMInfo$TOMSimilarities[[blockNo]] [, useIndivTOMSubset[set]]
        }
        if (useDiskCache)
        {
          # Order TOM (this may take a long time...)
          if (verbose > 3) printFlush(spaste(spaces, "......ordering TOM"));
          time = system.time({order1 = .qorder(tomDS)});
          if (verbose > 3) { printFlush("Time to order TOM:"); print(time); }
          # save the order
          orderFiles[set] = tempfile(pattern = spaste(".orderForSet", set), tmpdir = cacheDir);
          if (verbose > 3) printFlush(spaste(spaces, "......saving order and ordered TOM"));
          save(order1, file = orderFiles[set]);
          # Save ordered tomDS into chunks
          tomDS.ordered = tomDS[order1];
          sc = .saveChunks(tomDS.ordered, chunkSize, cacheBase, cacheDir = cacheDir);
          chunkFileNames[, set] = sc$files;
          chunkLengths = sc$chunkLengths;
        } else {
          setTomDS[, set] = tomDS[]
        }
      }
      if (useDiskCache)
      {
        # Step 2: Load chunks one by one and quantile normalize
        if (verbose > 2) printFlush(spaste(spaces, "....quantile normalizing chunks"));
        for (c in 1:nChunks)
        {
          if (verbose > 3) printFlush(spaste(spaces, "......QN for chunk ", c, " of ", nChunks));
          chunkData = matrix(NA, chunkLengths[c], nSets);
          for (set in 1:nSets)
            chunkData[, set] = .loadObject(chunkFileNames[c, set]);

          time = system.time({ chunk.norm = normalize.quantiles(chunkData, copy = FALSE);});
          if (verbose > 1) { printFlush("Time to QN chunk:"); print(time); }
          # Save quantile normalized chunks
          for (set in 1:nSets)
          {
            temp = chunk.norm[, set];
            save(temp, file = chunkFileNames[c, set]);
          }
        }

        if (verbose > 2) printFlush(spaste(spaces, "....putting together full QN'ed TOMs"));
        # Put together full TOMs
        for (set in 1:nSets)
        {
           load(orderFiles[set]);
           start = 1;
           for (c in 1:nChunks)
           {
             end = start + chunkLengths[c] - 1;
             tomDS[order1[start:end]] = .loadObject(chunkFileNames[c, set], size = chunkLengths[c]);
             start = start + chunkLengths[c];
           }
           if (saveCalibratedIndividualTOMs)
              save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo]);
           .saveChunks(tomDS, chunkSize, fileNames = chunkFileNames[, set]);
           unlink(orderFiles[set]);
        }
      } else {
        # If disk cache is not being used, simply call normalize.quantiles on the full set.
        setTomDS = normalize.quantiles(setTomDS);
        if (saveCalibratedIndividualTOMs) for (set in 1:nSets)
        {
           tomDS = .vector2dist(setTomDS[, set]);
           save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo]);
        }
      }
    } else stop("Unrecognized value of 'networkCalibration': ", networkCalibration);

    # Calculate consensus network
    if (verbose > 2)
      printFlush(paste(spaces, "....Calculating consensus network"));
    if (useDiskCache)
    {
      start = 1;
      for (chunk in 1:nChunks)
      {
        if (verbose > 3) printFlush(paste(spaces, "......working on chunk", chunk));
        end = start + chunkLengths[chunk] - 1;
        setChunks = array(0, dim = c(chunkLengths[chunk], nSets));
        for (set in 1:nSets)
        {
          load(file = chunkFileNames[chunk, set]);
          setChunks[, set] = temp;
          file.remove(chunkFileNames[chunk, set]);
        }
        if (useMean | consensusQuantile > 0)
        {
          consTomDS[start:end] = .consensusCalculation.base(
                     setChunks, useMean = useMean, setWeightMat = setWeightMat,
                     consensusQuantile = consensusQuantile)$consensus;
        } else {
          tmp = .consensusCalculation.base(setChunks, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile);
          consTomDS[start:end] = tmp$consensus;
          countIndex = as.numeric(names(tmp$originCount));
          originCount[countIndex] = originCount[countIndex] + tmp$originCount; 
          rm(tmp);
        } 
        start = end + 1;
      }
    } else {
      if (useMean | consensusQuantile > 0)
      {
         consTomDS[] = .consensusCalculation.base(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)$consensus;
      } else {
          tmp = .consensusCalculation.base(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile);
          consTomDS[] = tmp$consensus;
          countIndex = as.numeric(names(tmp$originCount));
          originCount[countIndex] = originCount[countIndex] + tmp$originCount; 
          rm(tmp);
      }
    }
    
    # Save the consensus TOM if requested

    if (saveConsensusTOMs)
    {
       TOMFiles[blockIndex] = .substituteTags(consensusTOMFilePattern, "%b", blockNo);
       if (TOMFiles[blockIndex]==consensusTOMFilePattern)
         stop(paste("File name for consensus TOM must contain the tag %b somewhere in the file name -\n",
                    "   - this tag will be replaced by the block number. "));
       save(consTomDS, file = TOMFiles[blockIndex]);
    }

    if (returnTOMs) consensusTOM.out[[blockIndex]] = consTomDS;

    gc();
  }

  if (!saveConsensusTOMs) TOMFiles = NULL;
  if (!returnTOMs)  consensusTOM.out = NULL;

  if (localIndividualTOMCalculation)
  {
     if (!individualTOMInfo$saveTOMs)
     {
        # individual TOMs are contained in the individualTOMInfo list; remove them.
        individualTOMInfo$TOMSimilarities = NULL;
     }
  } 
        
  list(consensusTOM = consensusTOM.out,
       TOMFiles = TOMFiles,
       saveConsensusTOMs = saveConsensusTOMs,

       individualTOMInfo = individualTOMInfo,
       useIndivTOMSubset = useIndivTOMSubset,
       goodSamplesAndGenes = gsg,
       nGGenes = nGoodGenes,
       nSets = nSets,

       saveCalibratedIndividualTOMs = saveCalibratedIndividualTOMs,
       calibratedIndividualTOMFileNames = calibratedIndividualTOMFileNames,
       networkCalibrationSamples = if (getNetworkCalibrationSamples) networkCalibrationSamples else NULL,

       consensusQuantile = consensusQuantile, 
       originCount = originCount

      )
}
