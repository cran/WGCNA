# Hierarchical consensus modules

hierarchicalConsensusModules = function(
    multiExpr, 
    multiWeights = NULL,

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

    iteratePruningAndMerging = FALSE,
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

  haveWeights = !is.null(multiWeights);
  .checkAndScaleMultiWeights(multiWeights, multiExpr, scaleByMax = FALSE);

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
         multiWeights = multiWeights,
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
  
  allLabels = mergedLabels = rep(0, nGenes);
  allLabelIndex = NULL;

  # Restrict data to goodSamples and goodGenes

  gsg = consensusTOMInfo$individualTOMInfo$blockInfo$goodSamplesAndGenes;

  if (!gsg$allOK)
  {
    multiExpr = mtd.subset(multiExpr, gsg$goodSamples, gsg$goodGenes);
    multiExpr.imputed = mtd.subset(multiExpr.imputed, gsg$goodSamples, gsg$goodGenes);
    if (haveWeights) multiWeights = mtd.subset(multiWeights, gsg$goodSamples, gsg$goodGenes);
  }

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
  goodGeneLabels = rep(0, nGGenes);
  # Here's where the analysis starts

  for (blockNo in 1:nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."));
    # Select block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]];
    nBlockGenes = length(block);

    selExpr = mtd.subset(multiExpr, , block);
    if (haveWeights) selWeights = mtd.subset(multiWeights, , block);
    errorOccurred = FALSE;
    consTomDS = BD.getData(consensusTOMInfo$consensusData, blockNo);
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
    } else {
      blockLabels[blockLabels>0] = blockLabels[blockLabels>0] + maxUsedLabel;
      maxUsedLabel = max(blockLabels);
      goodGeneLabels[block] = blockLabels;
    }
  }

  prune = try(pruneAndMergeConsensusModules(
     multiExpr = multiExpr,
     multiWeights = multiWeights,
     multiExpr.imputed = multiExpr.imputed,
     labels = goodGeneLabels,

     networkOptions = networkOptions,
     consensusTree = consensusTree,

     minModuleSize = minModuleSize,
     minCoreKME = minCoreKME, 
     minCoreKMESize = minCoreKMESize,
     minKMEtoStay = minKMEtoStay,

     # Module eigengene calculation options

     impute = impute,
     trapErrors = trapErrors,

     # Module merging options

     calibrateMergingSimilarities = calibrateMergingSimilarities,
     mergeCutHeight = mergeCutHeight,

     iterate = iteratePruningAndMerging,
     collectGarbage = collectGarbage,
     getDetails = TRUE,
     verbose = verbose, indent=indent), silent = TRUE);

  if (inherits(prune, "try-error"))
  {
    printFlush(paste(spaces, "'pruneAndMergeConsensusModules' failed with the following error message:\n",
                     spaces, prune, "\n", spaces, "--> returning unpruned module labels.")); 
    mergedLabels = goodGeneLabels;
  } else 
    mergedLabels = prune$labels;

  allLabels[gsg$goodGenes] = goodGeneLabels;

  MEs = try(multiSetMEs(multiExpr, universalColors = mergedLabels,
                            # trapErrors = TRUE, returnValidOnly = TRUE
                            ), silent = TRUE);
  if (class(MEs)=='try-error')
  {
    warning(paste('blockwiseConsensusModules: ME calculation failed with this message:\n     ',
          MEs, '---> returning empty module eigengenes'));
    allSampleMEs = NULL;
  } else {
    mergedLabels[gsg$goodGenes] = mergedLabels;
    index = lapply(gsg$goodSamples, function(gs)
    {
      out = rep(NA, length(gs));
      out[gs] = 1:sum(gs);
      out;
    });
    allSampleMEs = mtd.subset(MEs, index);
  }

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

#=====================================================================================================
#
# pruneAndMergeConsensusModules
#
#=====================================================================================================

pruneConsensusModules = function(
  multiExpr,
  multiWeights = NULL,
  multiExpr.imputed = NULL,
  MEs = NULL,
  labels,

  unassignedLabel = if (is.numeric(labels)) 0 else "grey",

  networkOptions,
  consensusTree,

  minModuleSize,
  minCoreKMESize = minModuleSize/3,
  minCoreKME = 0.5,
  minKMEtoStay = 0.2,

  # Module eigengene calculation options
  impute = TRUE,
  collectGarbage = FALSE,
  checkWeights = TRUE, 

  verbose = 1, indent=0)
{
  spaces = indentSpaces(indent);

  if (checkWeights) .checkAndScaleMultiWeights(multiWeights, multiExpr, scaleByMax = FALSE)

  oldLabels = labels;
  moduleIndex = sort(unique(labels));
  moduleIndex = moduleIndex[moduleIndex!=0];
  if (is.null(MEs))
  {
    # If multiExpr.imputed were not given, do not impute here, let moduleEigengenes do it since the imputation
    # should be faster there.
    if (is.null(multiExpr.imputed)) multiExpr.imputed = multiExpr
    MEs = multiSetMEs(multiExpr.imputed, universalColors = labels,
                  excludeGrey = TRUE, grey = unassignedLabel, impute = impute,
                  verbose = verbose-4, indent = indent + 3);
  } else {
    meSize = checkSets(MEs);
  }

  deleteModules = numeric(0);
  changedModules = numeric(0);

  # Check modules: make sure that of the genes present in the module, at least a minimum number
  # have a correlation with the eigengene higher than a given cutoff, and that all member genes have
  # the required minimum consensus KME

  if (verbose>0) 
    printFlush(paste(spaces, "..checking kME in consensus modules"));

  nSets = nSets(multiExpr);

  if (is.null(multiWeights)) multiWeights = .listRep(numeric(0), nSets) 

  KME = mtd.mapply(function(expr, weights, me, netOpt)
    {
     # printFlush("=============================================================");
     # print(netOpt$corOptions);
      haveWeights = length(dim(weights))==2;
      kme = do.call(netOpt$corFnc, 
                 c(list(x = expr, y = me, weights.x = weights), netOpt$corOptions));
      if (!grepl("signed", netOpt$networkType)) kme = abs(kme);
      kme;
    }, multiExpr, multiWeights, MEs, networkOptions, returnList = TRUE);

  consKME = simpleHierarchicalConsensusCalculation(KME, consensusTree);

  if (collectGarbage) gc();

  nMEs = checkSets(MEs)$nGenes;

  for (mod in 1:nMEs)
  {
    modGenes = (labels==moduleIndex[mod]);
    consKME1 = consKME[modGenes, mod];
    if (sum(consKME1>minCoreKME) < minCoreKMESize) 
    {
      labels[modGenes] = 0;
      deleteModules = union(deleteModules, mod);
      if (verbose>1) 
        printFlush(paste(spaces, "    ..deleting module ",moduleIndex[mod], 
                         ": of ", sum(modGenes), 
                   " total genes in the module only ",  sum(consKME1>minCoreKME), 
                   " have the requisite high correlation with the eigengene in all sets.", sep=""));
    } else if (sum(consKME1<minKMEtoStay)>0)
    {
      if (verbose > 1) 
        printFlush(paste(spaces, "    ..removing", sum(consKME1<minKMEtoStay),
                         "genes from module", moduleIndex[mod], "because their KME is too low."));
      labels[modGenes][consKME1 < minKMEtoStay] = 0;
      if (sum(labels[modGenes]>0) < minModuleSize) 
      {
        deleteModules = union(deleteModules, mod);
        labels[modGenes] = unassignedLabel;
        if (verbose>1) 
          printFlush(paste(spaces, "    ..deleting module ",moduleIndex[mod], 
                   ": not enough genes in the module after removal of low KME genes.", sep=""));
      } else {
        changedModules = union(changedModules, moduleIndex[mod]);
      }
    }
  }

  # Remove marked modules

  if (length(deleteModules) > 0)
  {
     for (set in 1:nSets) MEs[[set]]$data = MEs[[set]]$data[, -deleteModules, drop = FALSE];
     modGenes = labels %in% moduleIndex[deleteModules];
     labels[modGenes] = unassignedLabel;
     moduleIndex = moduleIndex[-deleteModules];
  }

  labels;
}

pruneAndMergeConsensusModules = function(
  multiExpr,
  multiWeights = NULL,
  multiExpr.imputed = NULL,
  labels,

  unassignedLabel = if (is.numeric(labels)) 0 else "grey",
  networkOptions,
  consensusTree,

  minModuleSize,
  minCoreKMESize = minModuleSize/3,
  minCoreKME = 0.5, 
  minKMEtoStay = 0.2,

  # Module eigengene calculation options

  impute = TRUE,
  trapErrors = FALSE,

  # Module merging options

  calibrateMergingSimilarities = FALSE,
  mergeCutHeight = 0.15,

  iterate = TRUE,
  collectGarbage = FALSE,
  getDetails = TRUE,
  verbose = 1, indent=0)
  
{
  spaces = indentSpaces(indent);
  if (is.null(multiExpr.imputed)) 
    multiExpr.imputed = mtd.apply(multiExpr, function(x) t(impute.knn(t(scale(x)))$data));

  .checkAndScaleMultiWeights(multiWeights, multiExpr, scaleByMax = FALSE);

  changed = TRUE;
  if (getDetails) details = list(originalLabels = labels);
  step = 0;
  while (changed)
  {
    step = step + 1;
    if (verbose > 0) printFlush(spaste(spaces, "step ", step));
    stepDetails = list();
    oldLabels = labels;
    if (verbose>1) printFlush(paste(spaces, "..pruning genes with low KME.."));
    labels = pruneConsensusModules(
       multiExpr,
       multiWeights = multiWeights,
       multiExpr.imputed = multiExpr.imputed,
       MEs = NULL,
       labels = labels,

       unassignedLabel = unassignedLabel,
       networkOptions = networkOptions,
       consensusTree = consensusTree,
       minModuleSize = minModuleSize,
       minCoreKME = minCoreKME,
       minCoreKMESize = minCoreKMESize,
       minKMEtoStay = minKMEtoStay,

       impute = impute,
       collectGarbage = collectGarbage,
       checkWeights = FALSE,
       verbose = verbose -2, indent = indent + 1)

    if (getDetails) stepDetails = c(stepDetails, list(prunedLabels = labels))
    #if (sum(labels>0)==0)
    #{
    #  if (verbose>1) 
    #    printFlush(paste(spaces, "  ..No significant modules left."));
    #  if (getDetails) details = c(details, stepDetails);
    #  break;
    #}

    # Merging needs to be called only if we're either in first iteration or if pruning acutally changed
    # modules.
    if (step==1 || any(labels!=oldLabels))
    {
      if (verbose>1) printFlush(paste(spaces, "..merging consensus modules that are too close.."));

      mergedMods = hierarchicalMergeCloseModules(multiExpr, labels = labels,
                                networkOptions = networkOptions, consensusTree = consensusTree,
                                calibrateMESimilarities = calibrateMergingSimilarities,
                                cutHeight = mergeCutHeight,
                                relabel = TRUE,
                                verbose = verbose-2, indent = indent + 1);
      if (getDetails) stepDetails = c(stepDetails, list(mergeInfo = mergedMods));
      labels = mergedMods$labels;
    }
    changed = !all(labels==oldLabels) & iterate;
    if (getDetails) details = c(details, list(stepDetails));
  }
  if (getDetails) 
  {
     names(details)[-1] = spaste("Iteration.", prependZeros(1:step));
     list(labels = labels, lastMergeInfo = mergedMods, details = details);
  } else labels;
}

