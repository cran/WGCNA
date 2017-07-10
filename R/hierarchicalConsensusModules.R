# Hierarchical consensus modules

hierarchicalConsensusModules = function(multiExpr, 

         # Optional: multiExpr wigth imputed missing data
         multiExpr.imputed = NULL,

         # Data checking options

         checkMissingData = TRUE,

         # Blocking options

         blocks = NULL, 
         maxBlockSize = 5000, 
         blockSizePenaltyPower = 5,
         nPreclusteringCenters = NULL,
         randomSeed = 12345,

         # ...or information needed to construct individual networks

         # Network construction options. This can be a single object of class NetworkOptions, or a multiData
         # structure of NetworkOptions objects, one per element of multiExpr.

         networkOptions,

         # Save individual TOMs?

         saveIndividualTOMs = TRUE,
         individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
         keepIndividualTOMs = FALSE,

         # Consensus calculation options

         consensusTree = NULL,  # if not given, the one in consensusTOMInfo will be used.

         # Return options
         saveConsensusTOM = TRUE,
         consensusTOMFilePattern = "consensusTOM-%a-Block%b.RData",

         # Keep the consensus? Note: I will not have an option to keep intermediate results here.
         keepConsensusTOM = saveConsensusTOM,

         # Internal handling of TOMs

         useDiskCache = NULL, chunkSize = NULL,
         cacheBase = ".blockConsModsCache",
         cacheDir = ".",

         # Alternative consensus TOM input from a previous calculation 

         consensusTOMInfo = NULL,

         # Basic tree cut options 

         deepSplit = 2, 
         detectCutHeight = 0.995, minModuleSize = 20,
         checkMinModuleSize = TRUE,

         # Advanced tree cut opyions

         maxCoreScatter = NULL, minGap = NULL,
         maxAbsCoreScatter = NULL, minAbsGap = NULL,
         minSplitHeight = NULL, minAbsSplitHeight = NULL,

         useBranchEigennodeDissim = FALSE,
         minBranchEigennodeDissim = mergeCutHeight,

         stabilityLabels = NULL,
         stabilityCriterion = c("Individual fraction", "Common fraction"),
         minStabilityDissim = NULL,

         pamStage = TRUE,  pamRespectsDendro = TRUE,

         # Gene joining and removal from a module, and module "significance" criteria
         # reassignThresholdPS = 1e-4, ## For now do not do gene reassignment - have to think more about how
         # to do it.

         minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
         minKMEtoStay = 0.2,

         # Module eigengene calculation options

         impute = TRUE,
         trapErrors = FALSE,

         # Module merging options

         calibrateMergingSimilarities = FALSE,
         mergeCutHeight = 0.15, 
                          
         # General options
         collectGarbage = TRUE,
         verbose = 2, indent = 0,
         ...)
{
  spaces = indentSpaces(indent);

  dataSize = checkSets(multiExpr);
  nSets = dataSize$nSets;
  nGenes = dataSize$nGenes;
  # nSamples = dataSize$nSamples;

  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed, add = FALSE);
    } 
    set.seed(randomSeed);
  }

  if (checkMinModuleSize & (minModuleSize > nGenes/5))
  {
    minModuleSize = nGenes/5;
    warning("blockwiseConsensusModules: minModuleSize appeared too large and was lowered to", 
            minModuleSize,
            ". If you insist on keeping the original minModuleSize, set checkMinModuleSize = FALSE.");
  }

  if (verbose>0) 
     printFlush(paste(spaces, "Calculating consensus modules and module eigengenes", 
                      "block-wise from all genes"));

  # prepare scaled and imputed multiExpr.
  multiExpr.scaled = mtd.apply(multiExpr, scale);
  hasMissing = unlist(multiData2list(mtd.apply(multiExpr, function(x) { any(is.na(x)) })));
  # Impute those that have missing data
  if (is.null(multiExpr.imputed)) {
     multiExpr.imputed = mtd.mapply(function(x, doImpute) 
                         { if (doImpute) t(impute.knn(t(x))$data) else x },
                                   multiExpr.scaled, hasMissing);
  } else {
    size.imp = checkSets(multiExpr.imputed);
    if (!isTRUE(all.equal(size.imp, dataSize)))
      stop("If given, 'multiExpr.imputed' must have the same dimensions in each set as 'multiExpr'.");
  }
  branchSplitFnc = NULL;
  minBranchDissimilarities = numeric(0);
  externalSplitFncNeedsDistance = logical(0);
  if (useBranchEigennodeDissim)
  {
    branchSplitFnc = "hierarchicalBranchEigengeneDissim";
    minBranchDissimilarities = minBranchEigennodeDissim;
    externalSplitFncNeedsDistance = FALSE;
  } 

  if (!is.null(stabilityLabels))
  {
    stabilityCriterion = match.arg(stabilityCriterion);
    branchSplitFnc = c(branchSplitFnc, 
           if (stabilityCriterion=="Individual fraction") 
               "branchSplitFromStabilityLabels.individualFraction" else "branchSplitFromStabilityLabels");
    minBranchDissimilarities = c(minBranchDissimilarities, minStabilityDissim);
    externalSplitFncNeedsDistance = c(externalSplitFncNeedsDistance, FALSE);
  }

  otherArgs = list(...);

  getDetails = FALSE;
  if ("getDetails" %in% names(otherArgs)) getDetails = otherArgs$getDetails;


  if (is.null(consensusTOMInfo))
  {
    if (is.null(consensusTree))
      stop("Either 'consensusTOMInfo' or 'consensusTree' must be given.");

    consensusTOMInfo = hierarchicalConsensusTOM(
         multiExpr = multiExpr,
         checkMissingData = checkMissingData,
         blocks = blocks,
         maxBlockSize = maxBlockSize,
         blockSizePenaltyPower = blockSizePenaltyPower,
         nPreclusteringCenters = nPreclusteringCenters,
         randomSeed = NULL,

         networkOptions = networkOptions,

         keepIndividualTOMs = keepIndividualTOMs,
         individualTOMFileNames = individualTOMFileNames,

         consensusTree = consensusTree,

         saveCalibratedIndividualTOMs = FALSE,
         getCalibrationSamples = FALSE,

         # Return options
         saveConsensusTOM = saveConsensusTOM,
         consensusTOMFilePattern = consensusTOMFilePattern,

         keepIntermediateResults = FALSE,

         # Internal handling of TOMs
         useDiskCache = useDiskCache,
         chunkSize = chunkSize,
         cacheBase = cacheBase,
         cacheDir = cacheDir,
         collectGarbage = collectGarbage,
         verbose = verbose, indent = indent);

     removeConsensusTOMOnExit = !keepConsensusTOM;
  } else {
    # Basic checks on consensusTOMInfo
    .checkComponents(consensusTOMInfo, c("individualTOMInfo", "consensusData", "consensusTree"));

    if (length(consensusTOMInfo$individualTOMInfo$blockInfo$blocks)!=nGenes)
      stop("Inconsistent number of genes in 'consensusTOMInfo$individualTOMInfo$blockInfo$blocks'.");

    if (!is.null(consensusTree) && !isTRUE(all.equal(consensusTree, consensusTOMInfo$consensusTree)))
       warning(immediate. = TRUE,
              "hierarchicalConsensusModules: given 'consensusTree' is different\n",
              "from the 'consensusTree' component in 'consensusTOMInfo'. \n",
              "This is normally undesirable and may\n",
              "indicate a mistake in the function call.");

    if (is.null(consensusTree)) consensusTree = consensusTOMInfo$consensusTree;

    removeConsensusTOMOnExit = FALSE;
    networkOptions = consensusTOMInfo$individualTOMInfo$networkOptions;
  }
  
  allLabels = rep(0, nGenes);
  allLabelIndex = NULL;

  # Restrict data to goodSamples and goodGenes

  gsg = consensusTOMInfo$individualTOMInfo$blockInfo$goodSamplesAndGenes;

  if (!gsg$allOK)
    multiExpr = mtd.subset(multiExpr, gsg$goodSamples, gsg$goodGenes);

  nGGenes = sum(gsg$goodGenes);
  nGSamples = sapply(gsg$goodSamples, sum);

  blocks = consensusTOMInfo$individualTOMInfo$blockInfo$blocks;
  gBlocks = consensusTOMInfo$individualTOMInfo$blockInfo$gBlocks;

  blockLevels = sort(unique(gBlocks));
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels);

  # reassignThreshold = reassignThresholdPS^nSets;

  consMEs = vector(mode = "list", length = nSets);
  dendros = list();

  cutreeLabels = list();
  maxUsedLabel = 0;
  # Here's where the analysis starts

  for (blockNo in 1:nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."));
    # Select block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]];
    nBlockGenes = length(block);

    selExpr = mtd.subset(multiExpr, , block);
    errorOccurred = FALSE;
    consTomDS = BD.getData(consensusTOMInfo$consensusData, blockNo);
    # Temporary "cast" so fastcluster::hclust doesn't complain about non-integer size.
    # attr(consTomDS, "Size") = as.integer(attr(consTomDS, "Size")); ## This should not be needed now.

    consTomDS = 1-consTomDS;
    
    if (collectGarbage) gc();

    if (verbose>2) printFlush(paste(spaces, "....clustering and detecting modules.."));
    errorOccured = FALSE;
    dendros[[blockNo]] = fastcluster::hclust(consTomDS, method = "average");
    if (verbose > 8)
    {
      if (interactive())
        plot(dendros[[blockNo]], labels = FALSE, main = paste("Block", blockNo));
    }

    externalSplitOptions = list();
    e.index = 1;
    if (useBranchEigennodeDissim)
    {
      externalSplitOptions[[e.index]] = list(multiExpr = mtd.subset(multiExpr.imputed,, block),
               networkOptions = networkOptions,
               consensusTree = consensusTree);
      e.index = e.index +1;
    }
    if (!is.null(stabilityLabels))
    {
      externalSplitOptions[[e.index]] = list(stabilityLabels = stabilityLabels);
      e.index = e.index + 1;
    }

    #blockLabels = try(cutreeDynamic(dendro = dendros[[blockNo]], 
    blockLabels = cutreeDynamic(dendro = dendros[[blockNo]], 
                    distM = as.matrix(consTomDS), 
                    deepSplit = deepSplit,
                    cutHeight = detectCutHeight, minClusterSize = minModuleSize, 
                    method ="hybrid", 
                    maxCoreScatter = maxCoreScatter, minGap = minGap,
                    maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                    minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight,

                    externalBranchSplitFnc = branchSplitFnc, 
                    minExternalSplit = minBranchDissimilarities,
                    externalSplitOptions = externalSplitOptions,
                    externalSplitFncNeedsDistance = externalSplitFncNeedsDistance,
                    assumeSimpleExternalSpecification = FALSE,

                    pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                    verbose = verbose-3, indent = indent + 2)
                    #verbose = verbose-3, indent = indent + 2), silent = TRUE);
    if (verbose > 8)
    {
      print(table(blockLabels));
      if (interactive())
        plotDendroAndColors(dendros[[blockNo]], labels2colors(blockLabels), dendroLabels = FALSE, 
           main = paste("Block", blockNo));
    }
    if (getDetails)
    {
       cutreeLabels[[blockNo]] = blockLabels;
    }
    if (class(blockLabels)=='try-error')
    {
      (if (verbose>0) printFlush else warning)
           (paste(spaces, "blockwiseConsensusModules: cutreeDynamic failed:\n    ", spaces, 
                  blockLabels, "\n", spaces, "    Error occured in block", blockNo, "\n",
                  spaces, "   Continuing with next block. "));
      next;
    } else {
      blockLabels[blockLabels>0] = blockLabels[blockLabels>0] + maxUsedLabel;
      maxUsedLabel = max(blockLabels);
    }
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1) 
      {
          printFlush(paste(spaces, "No modules detected in block", blockNo,
                           "--> continuing with next block."))
      }
      next;
    }

    # Calculate eigengenes for this batch

    if (verbose>2) printFlush(paste(spaces, "....calculating eigengenes.."));
    blockAssigned = c(1:nBlockGenes)[blockLabels!=0];
    blockLabelIndex = sort(unique(blockLabels[blockAssigned]));
    blockConsMEs = try(multiSetMEs(selExpr, universalColors = blockLabels,
                                   excludeGrey = TRUE, grey = 0, impute = impute,
                                   # trapErrors = TRUE, returnValidOnly = TRUE, 
                                   verbose = verbose-4, indent = indent + 3), silent = TRUE);
    if (class(blockConsMEs)=='try-error')
    {
      if (verbose>0)
      {
        printFlush(paste(spaces, "*** multiSetMEs failed with the message:"));
        printFlush(paste(spaces, "     ", blockConsMEs));
        printFlush(paste(spaces, "*** --> Ending module detection here"));
      } else warning(paste("blocwiseConsensusModules: multiSetMEs failed with the message: \n",
               "      ", blockConsMEs, "\n--> continuing with next block."));
      next;
    }

    deleteModules = NULL;
    changedModules = NULL;

    if (collectGarbage) gc();

    # Check modules: make sure that of the genes present in the module, at least a minimum number
    # have a correlation with the eigengene higher than a given cutoff, and that all member genes have
    # the required minimum consensus KME

    if (verbose>2) 
      printFlush(paste(spaces, "....checking consensus modules for statistical meaningfulness.."));

    KME = mtd.mapply(function(expr, me, netOpt)
      {
       # printFlush("=============================================================");
       # print(netOpt$corOptions);
        kme = do.call(netOpt$corFnc, c(list(x = expr, y = me), netOpt$corOptions));
        if (!grepl("signed", netOpt$networkType)) kme = abs(kme);
        kme;
      }, selExpr, blockConsMEs, networkOptions, returnList = TRUE);
    consKME = simpleHierarchicalConsensusCalculation(KME, consensusTree);

    for (mod in 1:ncol(blockConsMEs[[1]]$data))
    {
      modGenes = (blockLabels==blockLabelIndex[mod]);
      consKME1 = consKME[modGenes, mod];
      if (sum(consKME1>minCoreKME) < minCoreKMESize) 
      {
        blockLabels[modGenes] = 0;
        deleteModules = union(deleteModules, mod);
        if (verbose>3) 
          printFlush(paste(spaces, "    ..deleting module ",blockLabelIndex[mod], 
                           ": of ", sum(modGenes), 
                     " total genes in the module only ",  sum(consKME1>minCoreKME), 
                     " have the requisite high correlation with the eigengene in all sets.", sep=""));
      } else if (sum(consKME1<minKMEtoStay)>0)
      {
        if (verbose > 3) 
          printFlush(paste(spaces, "    ..removing", sum(consKME1<minKMEtoStay),
                           "genes from module", blockLabelIndex[mod], "because their KME is too low."));
        blockLabels[modGenes][consKME1 < minKMEtoStay] = 0;
        if (sum(blockLabels[modGenes]>0) < minModuleSize) 
        {
          deleteModules = union(deleteModules, mod);
          blockLabels[modGenes] = 0;
          if (verbose>3) 
            printFlush(paste(spaces, "    ..deleting module ",blockLabelIndex[mod], 
                     ": not enough genes in the module after removal of low KME genes.", sep=""));
        } else {
          changedModules = union(changedModules, blockLabelIndex[mod]);
        }
      }
    }

    # Remove marked modules

    if (!is.null(deleteModules)) 
    {
       for (set in 1:nSets) blockConsMEs[[set]]$data = blockConsMEs[[set]]$data[, -deleteModules];
       modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]));
       blockLabels[modGenes] = 0;
       modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]));
       allLabels[modAllGenes] = 0;
       blockLabelIndex = blockLabelIndex[-deleteModules];
    }

    # Check whether there's anything left
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1) 
      {
        printFlush(paste(spaces, "  ..No significant modules detected in block", blockNo))
        printFlush(paste(spaces, "  ..continuing with next block."));
      }
      next;
    }

    # Update module eigengenes

    for (set in 1:nSets) 
      if (is.null(dim(blockConsMEs[[set]]$data))) 
        dim(blockConsMEs[[set]]$data) = c(length(blockConsMEs[[set]]$data), 1);

    if (is.null(consMEs[[1]]))
    {
       for (set in 1:nSets) consMEs[[set]] = list(data = blockConsMEs[[set]]$data);
    } else for (set in 1:nSets)
       consMEs[[set]]$data = cbind(consMEs[[set]]$data, blockConsMEs[[set]]$data);

    # Update allLabels

    allLabelIndex = c(allLabelIndex, blockLabelIndex);
    allLabels[gsg$goodGenes][block[blockAssigned]] = blockLabels[blockAssigned];

    collectGarbage();
  
  }

  if (verbose>1) printFlush(paste(spaces, "..merging consensus modules that are too close.."));

  #print(table(allLabels));
  #print(is.numeric(allLabels))

  mergedLabels = rep(NA, nGenes);

  mergedMods = try(hierarchicalMergeCloseModules(multiExpr, allLabels[gsg$goodGenes],
                            networkOptions = networkOptions, consensusTree = consensusTree,
                            calibrateMESimilarities = calibrateMergingSimilarities,
                            cutHeight = mergeCutHeight,
                            relabel = TRUE,
                            verbose = verbose-2, indent = indent + 2), silent = TRUE); 
  if (class(mergedMods)=='try-error')
  {
    if (verbose>0) 
    {
      printFlush(paste(spaces, 'blockwiseConsensusModules: mergeCloseModule failed with this message:\n',
            spaces, '    ', mergedMods, spaces,
            '---> returning unmerged consensus modules'));
    } else warning(paste('blockwiseConsensusModules: mergeCloseModule failed with this message:\n     ',
                          mergedMods, '---> returning unmerged consensus modules'));
    MEs = try(multiSetMEs(multiExpr, universalColors = allLabels[gsg$goodGenes]
                          # trapErrors = TRUE, returnValidOnly = TRUE
                          ), silent = TRUE);
    if (class(MEs)=='try-error')
    {
      warning(paste('blockwiseConsensusModules: ME calculation failed with this message:\n     ',
            MEs, '---> returning empty module eigengenes'));
      allSampleMEs = NULL;
    } else {
      mergedLabels[gsg$goodGenes] = allLabels[gsg$goodGenes];
      allSampleMEs = vector(mode = "list", length = nSets);
      for (set in 1:nSets)
      {
        allSampleMEs[[set]] =
           list(data = as.data.frame(matrix(NA, nrow = nGSamples[set], ncol = ncol(MEs[[set]]$data))));
        allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = MEs[[set]]$data[,];
        names(allSampleMEs[[set]]$data) = names(MEs[[set]]$data);
      }
    }
  } else {
    mergedLabels[gsg$goodGenes] = mergedMods$labels;
    allSampleMEs = vector(mode = "list", length = nSets);
    for (set in 1:nSets)
    {
      allSampleMEs[[set]] = 
         list(data = as.data.frame(matrix(NA, nrow = nGSamples[set], 
                                          ncol = ncol(mergedMods$newMEs[[1]]$data))));
      allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = mergedMods$newMEs[[set]]$data[,];
      names(allSampleMEs[[set]]$data) = names(mergedMods$newMEs[[set]]$data);
    }
  }

  names(allSampleMEs) = names(multiExpr);

  if (removeConsensusTOMOnExit) 
  {
    BD.checkAndDeleteFiles(consensusTOMInfo$consensusData);
    consensusTOMInfo$consensusData = NULL;
  }

  list(labels = mergedLabels,
       unmergedLabels = allLabels,
       colors = labels2colors(mergedLabels),
       unmergedColors = labels2colors(allLabels),
       multiMEs = allSampleMEs, 
       dendrograms = dendros,
       consensusTOMInfo = consensusTOMInfo,
       blockInfo = consensusTOMInfo$individualTOMInfo$blockInfo,
       moduleIdentificationArguments = list(
         deepSplit = deepSplit,
         detectCutHeight = detectCutHeight,
         minModuleSize = minModuleSize,
         maxCoreScatter = maxCoreScatter,
         minGap = minGap,
         maxAbsCoreScatter = maxAbsCoreScatter,
         minAbsGap = minAbsGap,
         minSplitHeight = minAbsSplitHeight,
         useBranchEigennodeDissim = useBranchEigennodeDissim,
         minBranchEigennodeDissim = minBranchEigennodeDissim,
         minStabilityDissim = minStabilityDissim,
         pamStage = pamStage,
         pamRespectsDendro = pamRespectsDendro,
         minCoreKME = minCoreKME,
         minCoreKMESize = minCoreKMESize,
         minKMEtoStay = minKMEtoStay,
         calibrateMergingSimilarities = calibrateMergingSimilarities,
         mergeCutHeight = mergeCutHeight),
       details = if(getDetails) list(cutreeLabels = cutreeLabels) else NULL
      );
}
