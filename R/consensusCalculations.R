consensusCalculation = function(
      # a list or multiData structure of either numeric vectors (possibly arrays) or blockwiseAdj objects
      individualData,  

      consensusOptions,
      
      useBlocks = NULL,
      randomSeed = NULL,
      saveCalibratedIndividualData = FALSE,
      calibratedIndividualDataFilePattern = "calibratedIndividualData-%a-Set%s-Block%b.RData",

      # Return options: the data can be either saved or returned but not both.
      saveConsensusData = TRUE,
      consensusDataFileNames = "consensusData-%a-Block%b.RData",
      getCalibrationSamples= FALSE,

      # Internal handling of data
      useDiskCache = NULL, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      # Behaviour
      collectGarbage = FALSE,
      verbose = 1, indent = 0)
{
  nSets = length(individualData);

  if (! isMultiData(individualData))
     individualData = list2multiData(individualData);

  setNames = names(individualData);
  if (is.null(setNames)) setNames = rep("", nSets);

  blockwise = inherits(individualData[[1]]$data, "BlockwiseData");

  if (!blockwise)
  {
    blockDimnames = .mtd.checkDimConsistencyAndGetDimnames(individualData);
    blockLengths = length(individualData[[1]]$data);
    blockAttributes = attributes(individualData[[1]]$data);
    metaData = list();
  } else {
    blockLengths = BD.blockLengths(individualData[[1]]$data);
    blockAttributes = individualData[[1]]$data$attributes;
    metaData = BD.getMetaData(individualData[[1]]$data, blocks = 1);
  }
  nBlocks = length(blockLengths);

  spaces = indentSpaces(indent);

  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed);
    } 
    set.seed(randomSeed);
  }

  setWeights = consensusOptions$setWeights;
  if (is.null(setWeights)) setWeights = rep(1, nSets);

  if (length(setWeights)!=nSets)
    stop("Length of 'setWeights' must equal the number of sets.");

  if (any(!is.finite(setWeights)))
    stop("setWeights must all be finite.");

  setWeightMat = as.matrix(setWeights)/sum(setWeights)

  if (is.null(chunkSize)) chunkSize = as.integer(.largestBlockSize/(2*nSets))
  if (is.null(useDiskCache)) useDiskCache = .useDiskCache(individualData, chunkSize = chunkSize);

  # Initialize various variables

  if (getCalibrationSamples)
  {
    if (!consensusOptions$sampleForCalibration)
      stop(paste("Incompatible input options: calibrationSamples can only be returned", 
                 "if sampleForCalibration is TRUE."));
    calibrationSamples = list();
  }

  blockLevels = 1:nBlocks;
  if (is.null(useBlocks)) useBlocks = blockLevels;
  useBlockIndex = match(useBlocks, blockLevels);

  if (!all(useBlocks %in% blockLevels))
    stop("All entries of 'useBlocks' must be valid block levels.");

  if (any(duplicated(useBlocks)))
    stop("Entries of 'useBlocks' must be unique.");

  nUseBlocks = length(useBlocks);
  if (nUseBlocks==0)
    stop("'useBlocks' cannot be non-NULL and empty at the same time.");

  consensus.out = list();

  consensusFiles = rep("", nUseBlocks);
  originCount = rep(0, nSets);

  calibratedIndividualDataFileNames = NULL;
  if (saveCalibratedIndividualData)
  {
    calibratedIndividualDataFileNames = matrix("", nSets, nBlocks);
    for (set in 1:nSets) for (b in 1:nBlocks)
      calibratedIndividualDataFileNames[set, b] = 
                .processFileName(calibratedIndividualDataFilePattern, setNumber = set, setNames = setNames,
                           blockNumber = b, analysisName = consensusOptions$analysisName);
  }
  if (collectGarbage) gc();

  calibratedIndividualData.saved = vector(mode = "list", length = nSets);
  consensusData = NULL;
  dataFiles = character(nUseBlocks);

  # Here's where the analysis starts
  for (blockIndex in 1:nUseBlocks)
  {
    block = useBlockIndex[blockIndex];

    if (verbose>1) printFlush(spaste(spaces, "..Working on block ", block, "."));

    scaleQuant = rep(1, nSets);
    scalePowers = rep(1, nSets);

    useDiskCache1 = useDiskCache && nSets > 1;  ### No need to use disk cache when there is only 1 set.
    # Set up file names or memory space to hold the set Data
    if (useDiskCache1)
    {
      nChunks = ceiling(blockLengths[block]/chunkSize);
      chunkFileNames = array("", dim = c(nChunks, nSets));
      on.exit(.checkAndDelete(chunkFileNames), add = TRUE);
    } else nChunks = 1;

    if (nChunks==1) useDiskCache1 = FALSE;
    if (!useDiskCache1)
    {
      # Note: setTomDS will contained the scaled set Data matrices.
      calibratedData = array(0, dim = c(blockLengths[block], nSets));
    } 

    # sample entry indices from the distance structure for Data scaling, if requested

    if (consensusOptions$calibration=="single quantile" && 
           consensusOptions$sampleForCalibration)
    {
      qx = min(consensusOptions$calibrationQuantile, 1-consensusOptions$calibrationQuantile);
      nScGenes = min(consensusOptions$sampleForCalibrationFactor * 1/qx, blockLengths[block]);
      scaleSample = sample(blockLengths[block], nScGenes);
      if (getCalibrationSamples)
        calibrationSamples[[blockIndex]] = list(sampleIndex = scaleSample,
                                            TOMSamples = matrix(NA, nScGenes, nSets));
    }
    if (consensusOptions$calibration %in% c("single quantile", "none"))
    {
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(spaste(spaces, "....Working on set ", set, " (", setNames[set], ")"))

        # We need to drop dimensions here but we may need them later. Keep that in mind.
        tomDS = as.numeric(BD.getData(individualData[[set]]$data, block, simplify = TRUE));
        
        if (consensusOptions$calibration=="single quantile")
        {
          # Scale Data so that calibrationQuantile agree in each set
          if (consensusOptions$sampleForCalibration)
          {
            if (getCalibrationSamples)
            { 
              calibrationSamples[[blockIndex]]$dataSamples[, set] = tomDS[scaleSample];
              scaleQuant[set] = quantile(calibrationSamples[[blockIndex]]$dataSamples[, set], 
                                         probs = consensusOptions$calibrationQuantile, type = 8);
            } else {
              scaleQuant[set] = quantile(tomDS[scaleSample], probs = consensusOptions$calibrationQuantile, 
                                         type = 8);
            }
          } else
            scaleQuant[set] = quantile(x = tomDS, probs = consensusOptions$calibrationQuantile, type = 8);
          if (set>1) 
          {
             scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
             tomDS = tomDS^scalePowers[set];
          }
          if (saveCalibratedIndividualData)
             calibratedIndividualData.saved[[set]] = 
                  addBlockToBlockwiseData(
                        calibratedIndividualData.saved[[set]],
                        .setAttrFromList(tomDS, blockAttributes[[blockIndex]]),
                        external = TRUE,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = calibratedIndividualDataFileNames[set, block])
        } 

        # Save the calculated Data either to disk in chunks or to memory.
      
        if (useDiskCache1)
        {
          if (verbose > 3) printFlush(paste(spaces, "......saving Data similarity to disk cache.."));
          sc = .saveChunks(tomDS, chunkSize, cacheBase, cacheDir = cacheDir);
          chunkFileNames[, set] = sc$files;
          chunkLengths = sc$chunkLengths;
        } else {
          calibratedData[, set] = tomDS;
        }
      }
      if (collectGarbage) gc();
    } else if (consensusOptions$calibration=="full quantile")
    {
      # Step 1: load each data set, get order, split Data into chunks according to order, and save.
      if (verbose>1) printFlush(spaste(spaces, "..working on quantile normalization"))
      if (useDiskCache1)
      {
        orderFiles = rep("", nSets);
        on.exit(.checkAndDelete(orderFiles),add = TRUE);
      }
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(spaste(spaces, "....Working on set ", set, " (", setNames[set], ")"))
        tomDS = as.numeric(BD.getData(individualData[[set]]$data, block, simplify = TRUE));

        if (useDiskCache1)
        {
          # Order Data (this may take a long time...)
          if (verbose > 3) printFlush(spaste(spaces, "......ordering Data"));
          time = system.time({order1 = .qorder(tomDS)});
          if (verbose > 1) { printFlush("Time to order Data:"); print(time); }
          # save the order
          orderFiles[set] = tempfile(pattern = spaste(".orderForSet", set), tmpdir = cacheDir);
          if (verbose > 3) printFlush(spaste(spaces, "......saving order and ordered Data"));
          save(order1, file = orderFiles[set]);
          # Save ordered tomDS into chunks
          tomDS.ordered = tomDS[order1];
          sc = .saveChunks(tomDS.ordered, chunkSize, cacheBase, cacheDir = cacheDir);
          chunkFileNames[, set] = sc$files;
          chunkLengths = sc$chunkLengths;
        } else {
           calibratedData[, set] = tomDS
        }
      }
      if (useDiskCache1)
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
        if (verbose > 2) printFlush(spaste(spaces, "....putting together full QN'ed Data"));
        # Put together full Data
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
           if (saveCalibratedIndividualData)
              calibratedIndividualData.saved[[set]] = addBlockToBlockwiseData(
                        calibratedIndividualData.saved[[set]],
                        .setAttrFromList(tomDS, blockAttributes[[blockIndex]]),
                        external = TRUE,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = calibratedIndividualDataFileNames[set, blockIndex]);
           .saveChunks(tomDS, chunkSize, fileNames = chunkFileNames[, set]);
           unlink(orderFiles[set]);
        }
      } else {
        # If disk cache is not being used, simply call normalize.quantiles on the full set.
        if (nSets > 1) calibratedData = normalize.quantiles(calibratedData);
        if (saveCalibratedIndividualData) for (set in 1:nSets)
        {
           calibratedIndividualData.saved[[set]] = addBlockToBlockwiseData(
                        calibratedIndividualData.saved[[set]],
                        .setAttrFromList(calibratedData[, set], blockAttributes[[blockIndex]]),
                        external = TRUE,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = calibratedIndividualDataFileNames[set, blockIndex]);

        }
      }
    } else stop("Unrecognized value of 'calibration' in consensusOptions: ", consensusOptions$calibration);

    # Calculate consensus 
    if (verbose > 2)
      printFlush(paste(spaces, "....Calculating consensus"));

    # create an empty consTomDS distance structure.
    consTomDS = numeric(blockLengths[block]);

    if (useDiskCache1)
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
        tmp = .consensusCalculation.base(setChunks, useMean = consensusOptions$useMean, 
                                         setWeightMat = setWeightMat,
                                         consensusQuantile = consensusOptions$consensusQuantile);
        consTomDS[start:end] = tmp$consensus;
        if (!is.null(tmp$originCount))
        {
          countIndex = as.numeric(names(tmp$originCount));
          originCount[countIndex] = originCount[countIndex] + tmp$originCount; 
        } 
        start = end + 1;
      }
    } else {
      tmp = .consensusCalculation.base(calibratedData, useMean = consensusOptions$useMean, 
                                       setWeightMat = setWeightMat,
                                       consensusQuantile = consensusOptions$consensusQuantile);
      consTomDS[] = tmp$consensus;
      if (!is.null(tmp$originCount))
      {
          countIndex = as.numeric(names(tmp$originCount));
          originCount[countIndex] = originCount[countIndex] + tmp$originCount; 
      }
    }
    
    # Save the consensus Data if requested
    if (saveConsensusData)
    {
       if (!grepl("%b", consensusDataFileNames))
         stop(paste("File name for consensus data must contain the tag %b somewhere in the file name -\n",
                    "   - this tag will be replaced by the block number. "));
       dataFiles[blockIndex] = .substituteTags(consensusDataFileNames, c("%b", "%a"), 
                                              c(block, consensusOptions$analysisName[1]));
    }
    consensusData = addBlockToBlockwiseData(
                        consensusData, 
                        .setAttrFromList(consTomDS, blockAttributes[[blockIndex]]),
                        external = saveConsensusData,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = if (saveConsensusData) dataFiles[blockIndex] else NULL)
    if (collectGarbage) gc();
  }

  list(
       #blockwiseData
       consensusData = consensusData,
       # int
       nSets = nSets,
       # Logical 
       saveCalibratedIndividualData = saveCalibratedIndividualData,
       # List of blockwise data of length nSets
       calibratedIndividualData = calibratedIndividualData.saved,
       # List with one component per block
       calibrationSamples = if (getCalibrationSamples) calibrationSamples else NULL,
       # Numeric vector with nSets components
       originCount = originCount,

       consensusOptions = consensusOptions

      )
}



#==========================================================================================================
#
# Hierarchical consensus calculation
#
#==========================================================================================================

#  hierarchical consensus tree: a list with the following components:
#    inputs: either an atomic character vector whose entries match names of individualData, or a list in
#      which each component can either be a single character string giving a name in individualDara, or
#      another hierarchical consensus tree.
#    consensusOptions: a list of class ConsensusOptions
#    analysisName: optional, analysis name used for naming files when saving to disk.

# Here individualData is a list or multiData in which every component is either a blockwiseData instance or
# a numeric object (matrix, vector etc). Function consensusCalculation handles both.

hierarchicalConsensusCalculation = function(
   individualData,
   
   consensusTree,

   level = 1,
   useBlocks = NULL,
   randomSeed = NULL,
   saveCalibratedIndividualData = FALSE,
   calibratedIndividualDataFilePattern = "calibratedIndividualData-%a-Set%s-Block%b.RData",

   # Return options: the data can be either saved or returned but not both.
   saveConsensusData = TRUE,
   consensusDataFileNames = "consensusData-%a-Block%b.RData",
   getCalibrationSamples= FALSE,

   # Return the intermediate results as well?
   keepIntermediateResults = FALSE,

   # Internal handling of data
   useDiskCache = NULL, chunkSize = NULL,
   cacheDir = ".",
   cacheBase = ".blockConsModsCache",

   # Behaviour
   collectGarbage = FALSE,
   verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);
  individualNames = names(individualData);
  if (is.null(individualNames)) 
    stop("'individualData' must be a named list.");

  if (!isMultiData(individualData)) individualData = list2multiData(individualData);

  if (!"inputs" %in% names(consensusTree))
    stop("'consensusTree' must contain component 'inputs'.");
  
  if (!"consensusOptions" %in% names(consensusTree))
    stop("'consensusTree' must contain component 'consensusOptions'.");

  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed);
    }
    set.seed(randomSeed);
  }

  # Set names for consensusTree$inputs so that the names are informative.

  if (is.null(names(consensusTree$inputs)))
  {
    names(consensusTree$inputs) = spaste("Level.", level, ".Input.", 
                                                  1:length(consensusTree$inputs));
    validInputNames = FALSE;
  } else 
    validInputNames = TRUE;

  isChar = sapply(consensusTree$inputs, is.character);
  names(consensusTree$inputs)[isChar] = consensusTree$inputs[isChar];
  if (!is.null(consensusTree$analysisName)) 
    consensusTree$consensusOptions$analysisName = consensusTree$analysisName;

  # Recursive step if necessary
  if (verbose > 0)
    printFlush(spaste(spaces, "------------------------------------------------------------------\n",
                      spaces, "   Working on ", consensusTree$consensusOptions$analysisName,
                      "\n", spaces, "------------------------------------------------------------------"));
  names(consensusTree$inputs) = make.unique(make.names(names(consensusTree$inputs)));
  inputs0 = mtd.mapply(function(inp1, name)
   {
     if (is.character(inp1))
     {
        if (!inp1 %in% names(individualData))
          stop("Element '", inp1, "' is not among names of 'individualData'.");
        inp1;
     } else {
        if ("analysisName" %in% names(inp1)) name1 = inp1$analysisName else name1 = name;
        inp1$consensusOptions$analysisName = name1;
        hierarchicalConsensusCalculation(individualData, inp1, 
                   useBlocks = useBlocks,
                   level = level + 1,
                   randomSeed = NULL,
                   saveCalibratedIndividualData = saveCalibratedIndividualData,
                   calibratedIndividualDataFilePattern =calibratedIndividualDataFilePattern,
                   saveConsensusData = saveConsensusData,
                   consensusDataFileNames = consensusDataFileNames,
                   getCalibrationSamples = getCalibrationSamples,
                   keepIntermediateResults = keepIntermediateResults,
                   useDiskCache = useDiskCache,
                   chunkSize = chunkSize,
                   cacheDir = cacheDir,
                   cacheBase = cacheBase,
                   collectGarbage = collectGarbage,
                   verbose = verbose -2, indent = indent + 2);
     }
  }, consensusTree$inputs, names(consensusTree$inputs));

  names(inputs0) = names(consensusTree$inputs)

  inputData = mtd.apply(inputs0, function(inp1)
  {
    if (is.character(inp1)) 
    {
       individualData[[inp1]]$data
    } else
       inp1$consensusData;
  });
    
  inputIsIntermediate = !sapply(consensusTree$inputs, is.character);

  # Need to check that all inputData have the same format. In particular, some could be plain numeric data and
  # some could be BlockwiseData.

  nInputs1 = length(inputData);
  isBlockwise = mtd.apply(inputData, inherits, "BlockwiseData", mdaSimplify = TRUE);
  if (any(!isBlockwise)) for (i in which(!isBlockwise))
     inputData[[i]]$data = newBlockwiseData(list(inputData[[i]]$data), external = FALSE)

  names(inputData) = names(consensusTree$inputs)

  # Calculate the consensus

  if (verbose > 0)
    printFlush(spaste(spaces, "..Final consensus calculation.."));
  consensus = consensusCalculation(
      individualData = inputData,
      consensusOptions = consensusTree$consensusOptions,
      randomSeed = NULL,
      saveCalibratedIndividualData = saveCalibratedIndividualData,
      calibratedIndividualDataFilePattern =calibratedIndividualDataFilePattern,
      saveConsensusData = saveConsensusData,
      consensusDataFileNames = consensusDataFileNames,
      getCalibrationSamples = getCalibrationSamples,
      useDiskCache = useDiskCache,
      chunkSize = chunkSize,
      cacheDir = cacheDir,
      cacheBase = cacheBase,
      collectGarbage = collectGarbage,
      verbose = verbose-1, indent = indent+1);

  if (saveConsensusData && !keepIntermediateResults && any(inputIsIntermediate))
    mtd.apply(inputData[inputIsIntermediate], BD.checkAndDeleteFiles);

  out = c(consensus, if (keepIntermediateResults) list(inputs = inputs0) else NULL); 

  out;
}


#==========================================================================================================
#
# Simple hierarchical consensus calculation from numeric data, with minimum checking and no calibration.
#
#==========================================================================================================

# Simpler version of consensus calculation, suitable for small data where calibration is not
# necessary.

simpleConsensusCalculation = function(
   # multiData or  list of numeric vectors 
   individualData,  
   consensusOptions,
   verbose = 1, indent = 0)
{
  nSets = length(individualData);

  if (isMultiData(individualData))
    individualData = multiData2list(individualData);

  if (consensusOptions$useMean)
  {
    setWeights = consensusOptions$setWeights;
    if (is.null(setWeights)) setWeights = rep(1, nSets);
    if (length(setWeights)!=nSets)
      stop("Length of 'setWeights' must equal the number of sets.");
  } else setWeights = NULL;

  .consensusCalculation.base.FromList(individualData, useMean = consensusOptions$useMean, 
                                   setWeights = setWeights,
                                   consensusQuantile = consensusOptions$consensusQuantile)$consensus;
}


# Simple hierarchical consensus

simpleHierarchicalConsensusCalculation = function(
   # multiData or  list of numeric vectors 
   individualData,
   consensusTree,
   level = 1)

{
  individualNames = names(individualData);
  if (is.null(individualNames)) 
    stop("'individualData' must be named.");

  if (is.null(names(consensusTree$inputs)))
    names(consensusTree$inputs) = spaste("Level.", level, ".Input.", 
                                                  1:length(consensusTree$inputs));

  if (isMultiData(individualData))
    individualData = multiData2list(individualData);

  isChar = sapply(consensusTree$inputs, is.character);
  names(consensusTree$inputs)[isChar] = consensusTree$inputs[isChar];

  # Recursive step if necessary

  names(consensusTree$inputs) = make.unique(make.names(names(consensusTree$inputs)));
  inputData = mapply(function(inp1, name)
   {
     if (is.character(inp1))
     {
        if (!inp1 %in% names(individualData))
          stop("Element '", inp1, "' is not among names of 'individualData'.");
        individualData[[inp1]];
     } else {
        if ("analysisName" %in% names(inp1)) name1 = inp1$analysisName else name1 = name;
        inp1$consensusOptions$analysisName = name1;
        simpleHierarchicalConsensusCalculation(individualData, inp1, 
                   level = level + 1)
     }
  }, consensusTree$inputs, names(consensusTree$inputs), SIMPLIFY = FALSE);

  # Calculate the consensus

  simpleConsensusCalculation(
      individualData = inputData,
      consensusOptions = consensusTree$consensusOptions)
}

