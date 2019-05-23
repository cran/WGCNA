# Call blockwise modules several times with sampled data, collect the results
# -02: first entry in the result will be the base modules.

# Also: add functions for determining module stability and whether modules should be merged. 

#===================================================================================================
#
# sampledBlockwiseModules
#
#===================================================================================================

sampledBlockwiseModules = function(
  datExpr, 
  nRuns, 
  startRunIndex = 1,
  endRunIndex = startRunIndex + nRuns -1,
  replace = FALSE, 
  fraction = if (replace) 1.0 else 0.63,
  randomSeed = 12345,
  checkSoftPower = TRUE,
  nPowerCheckSamples = 2000,
  skipUnsampledCalculation = FALSE,
  corType = "pearson",
  power = 6,
  networkType = "unsigned",
  saveTOMs = FALSE,
  saveTOMFileBase = "TOM",
  ...,
  verbose = 2, indent = 0)

{

  spaces = indentSpaces(indent);

  result = list();
  runTOMFileBase = saveTOMFileBase;
  nSamples = nrow(datExpr);
  nGenes = ncol(datExpr);

  corTypeI = pmatch(corType, .corTypes);
  if (is.na(corTypeI))
    stop(paste("Invalid 'corType'. Recognized values are", paste(.corTypes, collapse = ", ")))

  corFnc = .corFnc[corTypeI];
  
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

  if (checkSoftPower)
  {
    if (verbose > 0) printFlush(paste(spaces, "...calculating reference mean adjacencies.."));
    useGenes = sample(nGenes, nPowerCheckSamples, replace = FALSE);
    adj = adjacency(datExpr[, useGenes], power = power, type = networkType,
                      corFnc = corFnc)
    refAdjMeans = mean(as.dist(adj));
  }
  
  for (run in startRunIndex:endRunIndex)
  {
    set.seed(randomSeed + 2*run + 1);
    if (verbose > 0) printFlush(paste(spaces, "...working on run", run, ".."));
    if (saveTOMs) 
      runTOMFileBase = paste(saveTOMFileBase, "-run-", run, sep = "");
  
    if (run > startRunIndex || skipUnsampledCalculation)
    {
      useSamples = sample(nSamples, as.integer(nSamples * fraction), replace = replace)
    } else
      useSamples = c(1:nSamples)

    if (verbose > 2)
    {
      printFlush(paste(spaces, "Using the following samples: "))
      print(useSamples);
    }
    samExpr = as.matrix(datExpr[useSamples, ]);
    samPowers = power;
    if (checkSoftPower)
    {
      if (verbose > 1) printFlush(paste(spaces, "  ...calculating mean adjacencies in sampled data.."));
      adj = adjacency(samExpr[, useGenes], power = power, type = networkType,
                      corFnc = corFnc)
      sampledAdjMeans = mean(as.dist(adj));
      samPowers = power * log( refAdjMeans) / log( sampledAdjMeans);
      if (!is.finite(samPowers)) samPowers = power;
    }

    mods = blockwiseModules(
      datExpr = samExpr, 
      randomSeed = NULL,
      power = samPowers,
      corType = corType,
      networkType = networkType,
      saveTOMs = saveTOMs,
      saveTOMFileBase = runTOMFileBase,
      ...,
      verbose = verbose-2, indent = indent+2)
   
    result[[run]] = list(mods = mods, samples = useSamples, powers = samPowers)
  }

  if (seedSaved) .Random.seed <<- savedSeed;

  result;
}

#===================================================================================================
#
# sampledHierarchicalConsensusModules
#
#===================================================================================================

sampledHierarchicalConsensusModules = function(
  multiExpr, 
  multiWeights = NULL,

  networkOptions,
  consensusTree,

  nRuns, 
  startRunIndex = 1,
  endRunIndex = startRunIndex + nRuns -1,
  replace = FALSE, 
  fraction = if (replace) 1.0 else 0.63,
  randomSeed = 12345,
  checkSoftPower = TRUE,
  nPowerCheckSamples = 2000,
  individualTOMFilePattern = "individualTOM-Run.%r-Set%s-Block.%b.RData",
  keepConsensusTOMs = FALSE,
  consensusTOMFilePattern = "consensusTOM-Run.%r-%a-Block.%b.RData",
  skipUnsampledCalculation = FALSE,
  ...,
  verbose = 2, indent = 0,
  saveRunningResults = TRUE,
  runningResultsFile = "results.tmp.RData")
{

  spaces = indentSpaces(indent);

  result = list();
  exprSize = checkSets(multiExpr);
  nSets = exprSize$nSets;
  nSamples = exprSize$nSamples;

  .checkAndScaleMultiWeights(multiWeights, multiExpr, scaleByMax = FALSE);

  if (inherits(networkOptions, "NetworkOptions"))
    networkOptions = list2multiData(replicate(nSets, networkOptions, simplify = FALSE));
  
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit({ .Random.seed <<- savedSeed }, add = FALSE);
    }
    set.seed(randomSeed);
  }

  powers = unlist(mtd.apply(networkOptions, getElement, "power"));

  if (nPowerCheckSamples > exprSize$nGenes) nPowerCheckSamples = exprSize$nGenes;

  if (checkSoftPower)
  {
    if (verbose > 0) printFlush(paste(spaces, "...calculating reference mean adjacencies.."));
    useGenes = sample(exprSize$nGenes, nPowerCheckSamples, replace = FALSE);
    refAdjMeans = rep(0, nSets);
    for (set in 1:nSets)
    {
      adj = adjacency(multiExpr[[set]]$data[, useGenes], 
                      weights = if (is.null(multiWeights)) NULL else multiWeights[[set]]$data[, useGenes],
                      power = networkOptions[[set]]$data$power, type = networkOptions[[set]]$data$networkType,
                      corFnc = networkOptions[[set]]$data$corFnc,
                      corOptions = networkOptions[[set]]$data$corOptions)
      refAdjMeans[set] = mean(as.dist(adj));
    }
  }
  
  for (run in startRunIndex:endRunIndex)
  {
    runTOMFileBase = .substituteTags(consensusTOMFilePattern, "%r", run);
    individualTOMFiles1 = .substituteTags(individualTOMFilePattern, "%r", run);
    set.seed(randomSeed + 2*run + 1);
    
    if (verbose > 0) printFlush(paste(spaces, "Working on run", run, ".."));

    useSamples = list()
    for (set in 1:nSets)
    {
      if (run > startRunIndex-skipUnsampledCalculation)
      {
         printFlush("This run will be on sampled data.");
         useSamples[[set]] = sample(nSamples[set], as.integer(nSamples[set] * fraction), replace = replace)
      } else 
         useSamples[[set]] = c(1:nSamples[set]);
    }
    samExpr = mtd.subset(multiExpr, useSamples);
    if (!is.null(multiWeights)) 
    {
      samWeights = mtd.subset(multiWeights, useSamples) 
    } else {
      samWeights = NULL;
    }
    samPowers = powers;
    if (checkSoftPower)
    {
      if (verbose > 1) printFlush(paste(spaces, "  ...calculating mean adjacencies in sampled data.."));
      sampledAdjMeans = rep(0, nSets);
      for (set in 1:nSets)
      {
        adj = adjacency(samExpr[[set]]$data[, useGenes], 
                        weights = if (is.null(multiWeights)) NULL else samWeights[[set]]$data[, useGenes],
                        power = networkOptions[[set]]$data$power, type = networkOptions[[set]]$data$networkType,
                        corFnc = networkOptions[[set]]$data$corFnc,
                        corOptions = networkOptions[[set]]$data$corOptions)
        sampledAdjMeans[set] = mean(as.dist(adj));
      }
      samPowers = powers * log( refAdjMeans) / log( sampledAdjMeans);
      samPowers[!is.finite(samPowers)] = powers[!is.finite(samPowers)];
    }

    networkOptions1 = mtd.mapply(function(x, power) { x$power = power; x; },
                       networkOptions, samPowers);

    collectGarbage();

    mods = hierarchicalConsensusModules(
      multiExpr = samExpr, 
      multiWeights = samWeights,
      randomSeed = NULL,

      networkOptions = networkOptions1,
      consensusTree = consensusTree,

      consensusTOMFilePattern = runTOMFileBase,
      individualTOMFileNames = individualTOMFiles1,
                                           
      keepIndividualTOMs = FALSE,
      keepConsensusTOM = keepConsensusTOMs,
      ...,
      verbose = verbose-2, indent = indent+2)

    result[[run - startRunIndex + 1]] = list(mods = mods, samples = useSamples, powers = samPowers)

    if (saveRunningResults) save(result, file = runningResultsFile);

    print(lapply(mods, object.size));
    print(gc());

  }

  result;
}


  
