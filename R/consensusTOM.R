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
   rm(temp); collectGarbage();
   list(files = fileNames, chunkLengths = chunkLengths);
}

.loadObject = function(file, name = NULL, size = NULL)
{
  x = load(file);
  if (!is.null(name) && (x!=name))
    stop("File ", file, " does not contain object '", name, "'.")

  obj = get(x);
  if (!is.null(size) && (length(obj)!=size))
    stop("Object '", name, "' does not have the correct length.");
  obj;
}

.qorder = function(data)
{
  data = as.numeric(data);
  .Call("qorder", data)
}

# Actual consensus calculation distilled into one function. setTomMat is assumed to have sets in columns
# and gene pairs in rows. setWeightMat should be a matrix of dimensions (nSets, 1) and be normalized to sum=1.

.consensusCalculation = function(setTomMat, useMean, setWeightMat, consensusQuantile)
{
  if (useMean)
  {
     if (any(is.na(setTomMat)))
     {
       finiteMat = 1-is.na(setTomMat);
       setTomMat[is.na(setTomMat)] = 0;
       out = setTomMat %*% setWeightMat / finiteMat%*%setWeightMat;
     } else {
       out = setTomMat %*% setWeightMat;
     }
     out.list = list(consensus = out);
  } else if (consensusQuantile == 0) 
  {
      min = rep(0, nrow(setTomMat));
      which = rep(0, nrow(setTomMat));
      whichmin = .C("minWhichMin_row", as.double(setTomMat),
                    as.integer(nrow(setTomMat)), as.integer(ncol(setTomMat)),
                    as.double(min), as.double(which), PACKAGE = "WGCNA");
      min = whichmin[[4]];
      which = whichmin[[5]] + 1;
      rm(whichmin);
      out.list = list(consensus = min, originCount = table(as.integer(which)));
  } else {
     out.list = list(consensus = rowQuantileC(setTomMat, p = consensusQuantile));
  }
  out.list;
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
      consensusTOMFileNames = "consensusTOM-Block%b.RData",
      returnTOMs = FALSE,

      # Internal handling of TOMs
      useDiskCache = TRUE, chunkSize = NULL,
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
       seedSaved = TRUE;
       savedSeed = .Random.seed
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
      calibratedIndividualTOMFileNames[set, b] = .processFileName(calibratedIndividualTOMFilePattern, set, 
                                                              individualTOMInfo$setNames, b);
  }
  collectGarbage();

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
        rm(tomDS); collectGarbage();
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
          if (verbose > 1) { printFlush("Time to order TOM:"); print(time); }
          # save the order
          orderFiles[set] = tempfile(pattern = spaste("orderForSet", set), tmpdir = cacheDir);
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
          consTomDS[start:end] = .consensusCalculation(setChunks, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)$consensus;
        } else {
          tmp = .consensusCalculation(setChunks, useMean = useMean, setWeightMat = setWeightMat,
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
         consTomDS[] = .consensusCalculation(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)$consensus;
      } else {
          tmp = .consensusCalculation(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
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
       TOMFiles[blockIndex] = .substituteTags(consensusTOMFileNames, "%b", blockNo);
       if (TOMFiles[blockIndex]==consensusTOMFileNames)
         stop(paste("File name for consensus TOM must contain the tag %b somewhere in the file name -\n",
                    "   - this tag will be replaced by the block number. "));
       save(consTomDS, file = TOMFiles[blockIndex]);
    }

    if (returnTOMs) consensusTOM.out[[blockIndex]] = consTomDS;

    collectGarbage();
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
        

  if (seedSaved) .Random.seed <<- savedSeed;

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

