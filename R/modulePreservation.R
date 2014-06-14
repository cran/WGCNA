# module preservation for networks specified by adjacency matrices, not expression data.
# this version can handle both expression and adjacency matrices.

# 
# Changelog:

# 2010/08/04: 
#   . Adding clusterCoeff and MAR to density preservation statistics
#   . Adding cor.clusterCoeff and cor.MAR to connectivity preservation statistics
#   . Adding silhuette width to separability statistics

# Adding p-values to output. For each Z score can add a corresponding p-value, bonferoni corrected p-value,
# and q value (local FDR)


#=====================================================================================================
#
# p-value functions
#
#=====================================================================================================
# Assumes Z in the form of a list, Z[[ref]][[test]] is a matrix of Z scores except the first column is
# assumed to be moduleSize. 
# Returned value is log10 of the p-value

.pValueFromZ = function(Z, bonf = FALSE, summaryCols = NULL, summaryInd = NULL)
{
  p = Z # This carries over all necessary names and missing values when ref==test
  nRef = length(Z);
  for (ref in 1:nRef)
  {
    nTest = length(Z[[ref]])
    for (test in 1:nTest)
    {
      zz = Z[[ref]][[test]];
      if (length(zz) > 1)
      {
        names = colnames(zz);
        # Try to be a bit intelligent about whether to ignore the first column
        if (names[1]=="moduleSize") ignoreFirst = TRUE else ignoreFirst = FALSE;
        ncol = ncol(zz);
        range = c( (1+as.numeric(ignoreFirst)):ncol);
        p[[ref]][[test]][, range] = pnorm(as.matrix(zz[, range]), lower.tail = FALSE, log.p = TRUE)/log(10);
        Znames = names[range]
        if (bonf)
        {
           pnames = sub("Z", "log.p.Bonf", Znames, fixed= TRUE);
           p[[ref]][[test]][, range] = p[[ref]][[test]][, range] + log10( nrow(p[[ref]][[test]]));
           biggerThan1 = p[[ref]][[test]][, range] > 0;
           p[[ref]][[test]][, range][biggerThan1] = 0; # Remember that the p-values are stored in log form
        } else
           pnames = sub("Z", "log.p", Znames, fixed= TRUE);
        if (!is.null(summaryCols)) for (c in 1:length(summaryCols))
        {
          medians = apply(p[[ref]][[test]][, summaryInd[[c]], drop = FALSE], 1, median, na.rm = TRUE)
          p[[ref]][[test]][,summaryCols[c]] = medians
        }
        colnames(p[[ref]][[test]])[range] = pnames;
      } 
    }
  }
  p;
}

.qValueFromP = function(p, summaryCols = NULL, summaryInd = NULL)
{
  q = p # This carries over all necessary names and missing values when ref==test
  nRef = length(p);
  for (ref in 1:nRef)
  {
    nTest = length(p[[ref]])
    for (test in 1:nTest)
    {
      pp = p[[ref]][[test]];
      if (length(pp) > 1)
      {
        names = colnames(pp);
        # Try to be a bit intelligent about whether to ignore the first column
        if (names[1]=="moduleSize")
        {
           ignoreFirst = TRUE 
        } else
           ignoreFirst = FALSE;
        ncol = ncol(pp);
        nrow = nrow(pp);
        range = c( (1+as.numeric(ignoreFirst)):ncol);
        for (col in range)
        {
          xx = try(qvalue(10^pp[,col]), silent = TRUE)
          if (class(xx)=="try-error" || length(xx)==1)
          {
            q[[ref]][[test]][, col] = rep(NA, nrow);
            printFlush(paste("Warning in modulePreservation: qvalue calculation failed for",
                             "column", col, "for reference set", ref, "and test set", test))
          } else 
            q[[ref]][[test]][, col] = xx$qvalues;
        }
        colnames(q[[ref]][[test]])[range] = sub("log.p", "q", names[range], fixed= TRUE);
        if (!is.null(summaryCols)) for (c in 1:length(summaryCols))
        {
          xx = log10(as.matrix(q[[ref]][[test]][, summaryInd[[c]], drop = FALSE]));
          # Restrict all -logs > 1000 to 1000:
          xx[!is.na(xx) & !is.finite(xx)] = -1000;
          xx[!is.na(xx) & (xx < -1000)] = -1000;
          medians  = apply(xx, 1, median, na.rm = TRUE)
          q[[ref]][[test]][, summaryCols[c]] = 10^medians;
        }
      }
    }
  }
  q;
}


modulePreservation = function(
   multiData,
   multiColor,
   dataIsExpr = TRUE,
   networkType = "unsigned",
   corFnc = "cor",
   corOptions = "use = 'p'",
   referenceNetworks = 1,
   testNetworks = NULL,
   nPermutations = 100,
   includekMEallInSummary = FALSE,
   restrictSummaryForGeneralNetworks = TRUE,
   calculateQvalue = FALSE, 
   randomSeed = 12345, 
   maxGoldModuleSize = 1000, maxModuleSize = 1000, 
   quickCor = 1,
   ccTupletSize = 2,
   calculateCor.kIMall = FALSE, 
   calculateClusterCoeff = FALSE,
   useInterpolation = FALSE,
   checkData = TRUE,
   greyName = NULL,
   savePermutedStatistics = TRUE,
   loadPermutedStatistics = FALSE,
   permutedStatisticsFile = 
      if (useInterpolation) "permutedStats-intrModules.RData" else "permutedStats-actualModules.RData",
   plotInterpolation = TRUE,
   interpolationPlotFile = "modulePreservationInterpolationPlots.pdf",
   discardInvalidOutput = TRUE,
   parallelCalculation = FALSE,
   verbose = 1, indent = 0)
{
   sAF = options("stringsAsFactors")
   options(stringsAsFactors = FALSE);
   on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)

   if (!is.null(randomSeed))
   {
     if (exists(".Random.seed"))
     {
        savedSeed = .Random.seed
        on.exit({.Random.seed <<- savedSeed;}, TRUE);
     }
     set.seed(randomSeed);
   }

   spaces = indentSpaces(indent);

   nType = charmatch(networkType, .networkTypes);
   if (is.na(nType))
      stop(paste("Unrecognized networkType argument.",
           "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")));

   # Check that the multiData/multiAdj has correct structure

   nNets = length(multiData);
   nGenes = sapply(multiData, sapply, ncol);
   if (checkData)
   {
     if (dataIsExpr)
     {  
       .checkExpr(multiData, verbose, indent)
     } else {
       .checkAdj(multiData, verbose, indent)
     }
   }

   # Check for names; if there are none, create artificial labels.
   setNames = names(multiData);
   if (is.null(setNames))
   {
     setNames = paste("Set", c(1:nNets), sep="");
   } 

   # Check that referenceNetworks is valid

   referenceNetworks = as.numeric(referenceNetworks);
   if (any(is.na(referenceNetworks)))
      stop("All elements of referenceNetworks must be numeric and present.");

   if (any(referenceNetworks < 1) | any(referenceNetworks > nNets))
      stop("referenceNetworks contains elements outside of the allowed range. ");

   nRefNets = length(referenceNetworks);

   # Check testNetworks

   if (is.null(testNetworks))
   {
      testNetworks = list();
      for (ref in 1:nRefNets) testNetworks[[ref]] = c(1:nNets)[ -referenceNetworks[ref] ];
   }

   # Check that testNetworks was specified correctly
   if (!is.list(testNetworks) && nRefNets > 1) 
     stop("When there are more than 1 reference networks, 'testNetworks' must\n",
          "  be a list with one component per reference network.");

   if (!is.list(testNetworks)) testNetworks = list(testNetworks);

   if (length(testNetworks)!=nRefNets) 
     stop("Length of 'testNetworks' must be the same as length of 'referenceNetworks'.");

   for (ref in 1:nRefNets)
   {
     if (any(testNetworks[[ref]] < 1 || testNetworks[[ref]] > nNets))
       stop("Some entries of testNetworks[[", ref, "]] are out of range.");
   }
   
   # Check multiColor

   if (class(multiColor)!="list")
   {
       stop("multiColor does not appear to have the correct format.")
   }

   if (length(multiColor)!=nNets)
   {
     multiColor2=list()
     if (length(names(multiColor))!=length(multiColor))
       stop("Each entry of 'multiColor' must have a name.");
     color2expr = match(names(multiColor), setNames);
     if (any(is.na(color2expr)))
       stop("Entries of 'multiColor' must name-match entries in 'multiData'.");
     for(s in 1:nNets)
     {
        #multiData[[s]]$data = as.matrix(multiData[[s]]$data);
        loc = which(names(multiColor) %in% setNames[s])
        if (length(loc)==0)
        {
          multiColor2[[s]] = NA
        } else
          multiColor2[[s]]=multiColor[[loc]]
     }
     multiColor = multiColor2
     rm(multiColor2);
   }

   s = 1; while ( (s <= nNets) & !.cvPresent(multiColor[[s]])) s = s+1;
   if (s==nNets+1)
     stop("No valid color lists found.")

   if (is.null(greyName))
   {
     if (is.numeric(multiColor[[s]]))  # Caution: need a valid s here.
     {
       greyName = 0
       goldName = 0.1;
     } else {
       greyName = "grey"
       goldName = "gold"
     }
   } else  {
     if (is.numeric(greyName)) goldName = 0.1 else goldName = "gold"
   }

   if (verbose > 2) printFlush(paste("  ..unassigned 'module' name:", greyName, 
                                     "\n  ..all network sample 'module' name:", goldName));

   MEgrey = paste("ME", greyName, sep="");
   MEgold = paste("ME", goldName, sep="");

   for (s in 1:nNets)
   {
     if ( .cvPresent(multiColor[[s]]) & nGenes[s]!=length(multiColor[[s]]))
       stop(paste("Color vector for set", s, "does not have the correct number of entries."));
     if (is.factor(multiColor[[s]])) multiColor[[s]] = as.character(multiColor[[s]]);
   }

   keepGenes = list();
   nNAs = sapply(multiColor, .nNAColors)
   if (any(nNAs > 0))
   {
     if (verbose > 0) printFlush(paste(spaces, " ..removing genes with missing color..."))
     for (s in 1:nNets)
       if (.cvPresent(multiColor[[s]]))
       {
          keepGenes[[s]] = !is.na(multiColor[[s]]);
          if (dataIsExpr)
          {
             multiData[[s]]$data = multiData[[s]]$data[, keepGenes[[s]]];
          } else
             multiData[[s]]$data = multiData[[s]]$data[keepGenes[[s]], keepGenes[[s]]];
          multiColor[[s]] = multiColor[[s]][keepGenes[[s]]];
       }
   }

   # Check for set names; if there are none, create artificial labels.

   if (is.null(names(multiData)))
   {
     setNames = paste("Set_", c(1:nNets), sep="");
   } else {
     setNames = names(multiData);
   }

   # Check that multiData has valid colnames

   for (s in 1:nNets)
      if (is.null(colnames(multiData[[s]]$data))) 
         stop(paste("Matrix of data in set", names(multiData)[s], 
                    "has no colnames. Colnames are needed to match variables."));

   # For now we use numeric labels for permutations.
   permGoldName = 0.1;
   permGreyName = 0;

   collectGarbage();

   if (verbose > 0) printFlush(paste(spaces, " ..calculating observed preservation values"))
   observed = .modulePreservationInternal(multiData, multiColor, dataIsExpr = dataIsExpr,
                     calculatePermutation = FALSE, networkType = networkType,
                     referenceNetworks = referenceNetworks, 
                     testNetworks = testNetworks, 
                     densityOnly = useInterpolation,
                     maxGoldModuleSize = maxGoldModuleSize,
                     maxModuleSize = maxModuleSize, 
                     corFnc = corFnc, corOptions = corOptions, 
                     quickCor = quickCor,
                     # calculateQuality = calculateQuality,
                     ccTupletSize = ccTupletSize,
                     calculateCor.kIMall = calculateCor.kIMall,
                     calculateClusterCoeff = calculateClusterCoeff,
                     checkData = FALSE, greyName = greyName, 
                     verbose = verbose -3, indent = indent + 2);

   if (nPermutations==0) return(list(observed = observed))

   # Calculate preservation scores in permuted data.

   psLoaded = FALSE;
   if (loadPermutedStatistics)
   {
     cat(paste(spaces, "..attempting to load permutation statistics.."));
     x = try(load(file=permutedStatisticsFile), silent = TRUE);
     if (class(x)=="try-error")
     {
       printFlush(paste("failed. Error message returned by system:\n", x));
     } else {
       expectVars = c("regModuleSizes", "regStatNames", "regModuleNames", "fixStatNames", "fixModuleNames", 
                      "permutationsPresent", "permOut", "interpolationUsed");
       e2x = match(expectVars, x);
       if (any(is.na(e2x)) | length(x)!=length(expectVars))
       {
         printFlush(paste("the file does not contain (all) expected variables."));
       } else if (length(permOut) != length(referenceNetworks))
       {
         printFlush("the loaded permutation statistics have incorrect number of reference sets.");
       } else if (length(permOut[[1]]) != nNets)
       {
         printFlush("the loaded permutation statistics have incorrect number of test sets.");
       } else {
         psLoaded = TRUE;
         printFlush("success.");
       }
     }
     if (!psLoaded)
       printFlush(paste(spaces, "\n ..will recalculate permutations.", 
                                "Hit Ctrl-C (Esc in Windows) to stop the calculation."));
   }

   nRegStats = 20;
   nFixStats = 3;
   
   if (!psLoaded)
   {
      permOut=list()      
      regModuleSizes = list();
      regModuleNames = list();      permutationsPresent = matrix(FALSE, nNets, nRefNets);
      interpolationUsed = matrix(FALSE, nNets, nRefNets);
      if (verbose > 0) printFlush(paste(spaces, " ..calculating permutation Z scores"))
      for(iref in 1:nRefNets)   # Loop over reference networks
      {
         ref = referenceNetworks[iref]
         if (verbose > 0) printFlush(paste(spaces, "..Working with set", ref, "as reference set"));
         permOut[[iref]] = list()      
         regModuleSizes[[iref]] = list();
         regModuleNames[[iref]] = list();
         nRefMods = length(unique(multiColor[[ref]]));
         if (nRefMods==1)
         {
           printFlush(paste(spaces, "*+*+*+*+* Reference set contains a single module.\n", 
                            spaces, 
                            "A permutation analysis is not meaningful for a single module; skipping."))
           next;
         }
         for (tnet in 1:nNets) if (tnet %in% testNetworks[[iref]])
         {
            # Retain only genes that are shared between the reference and test networks
            if (verbose > 1) printFlush(paste(spaces, "....working with set", tnet, "as test set"));
            overlap=intersect(colnames(multiData[[ref]]$data),colnames(multiData[[tnet]]$data))
            loc1=match(overlap, colnames(multiData[[ref]]$data))
            loc2=match(overlap, colnames(multiData[[tnet]]$data))
            refName = paste("ref_", setNames[ref],sep="")
            colorRef = multiColor[[ref]][loc1]
            if (dataIsExpr)
            {
              datRef=multiData[[ref]]$data[ , loc1]
              datTest=multiData[[tnet]]$data[ , loc2]
            } else {
              datRef=multiData[[ref]]$data[loc1, loc1]
              datTest=multiData[[tnet]]$data[loc2, loc2]
            }
            testName=setNames[tnet]
            nRefGenes = ncol(datRef);
              
            #if(!is.na(multiColor[[tnet]][1])||length(multiColor[[tnet]])!=1)
            if (.cvPresent(multiColor[[tnet]]))
            {
               colorTest=multiColor[[tnet]][loc2]
            } else  {
               colorTest=NA 
            }
            name=paste(refName,"vs",testName,sep="")
            obsModSizes=list()
            nObsMods = rep(0, 2);
            tab = table(colorRef);
            nObsMods[1] = length(tab);
            obsModSizes[[1]]=tab[names(tab)!=greyName]
            if ( !useInterpolation | (nObsMods[1] <= 5) | (sum(obsModSizes[[1]]) < 1000))
            {
               # Do not use interpolation: simply use original colors
               permRefColors = colorRef;
               interpolationUsed[tnet, iref] = FALSE;
               nPermMods = nObsMods[1];
               if (useInterpolation && (verbose > 1)) 
                  printFlush(paste(spaces, "    FYI: interpolation will not be used for this comparison."))
            } else {
               obsModSizes[[1]][obsModSizes[[1]]<3]=3
               if(length(colorTest)>1)
               {
                  tab = table(colorTest);
                  obsModSizes[[2]]=tab[names(tab)!=greyName]
                  obsModSizes[[2]][obsModSizes[[2]]<3]=3
                  nObsMods[2] = length(tab);
               } else {
                  obsModSizes[[2]]=NA
               }
               # Note we only need permColors for the reference set.
               nPermMods = 10
               minNMods = 5;
               OMS=obsModSizes[[1]]
               logmin=log(min(OMS)); logmax= log(min(maxModuleSize,max( OMS)))
               ok = FALSE;
               skip = FALSE;
               while (!ok)
               {
                  if (logmin>= logmax) logmax=logmin+nPermMods/2;
                  permModSizes=as.integer(exp(seq(from=logmin, to=logmax, length=nPermMods )))
                  nNeededGenes = sum(permModSizes)
                  # Check that the data has enough genes to fit the module sizes:
                  if (nNeededGenes < nRefGenes)
                  {
                     ok = TRUE; skip = FALSE;
                  } else if (nPermMods > minNMods)
                  {
                     # Drop one module
                     nPermMods = nPermMods -1;
                     logStep = (logmax - logmin)/nPermMods;
                     # If the modules are spaced far enough apart, also decrease the max module size
                     if (logStep > log(2))
                     {
                        logmax = logmax - logStep;
                     } 
                  } else {
                     # It appears we don't have enough genes to form meaningful modules for interpolation. 
                     # Decrease logmin as well, but only to a certain degree. 
                     logminFloor = 3;
                     if (logmin > logminFloor)
                     {
                        logmin = min(logmin-1, logminFloor);
                     } else {
                        # For now we give up, but in the future may have to add a non-interpolation approach here
                        # as well.
                        printFlush(paste(spaces, 
                                    "*+*+*+*+*+ There are not enough genes and/or modules for",
                                    "reference set", ref, "and test set", tnet, ".\n",
                                    spaces, "Will skip this combination."))
                        ok = TRUE;
                        skip = TRUE;
                     }
                  }
               }
      
               if (skip) 
               {
                  # Instead of skipping, use the actual colors
                  permRefColors = colorRef;
                  interpolationUsed[tnet, iref] = FALSE;
                  nPermMods = nObsMods[1];
               } else {
                  # Create a base label sequence for the permuted reference data set
                  permRefColors = rep(c(1:nPermMods), permModSizes);
                  nGrey = nRefGenes - length(permRefColors);
                  permRefColors = c(permRefColors, rep(greyName, nGrey));
                  interpolationUsed[tnet, iref] = TRUE;
               }
            }
            
            permExpr=list()
            permExpr[[1]]=list(data=datRef)
            permExpr[[2]]=list(data=datTest)
            names(permExpr) = setNames[c(ref, tnet)];
            permOut[[iref]][[tnet]]=list(
                      regStats = array(NA, dim = c(nPermMods+2-(!interpolationUsed[tnet, iref]), 
                                                   nRegStats, nPermutations)),
                      fixStats = array(NA, dim = c(nObsMods[[1]], nFixStats, nPermutations)));
   
            # Perform actual permutations
            #oldRNG = NULL;
            permColors = list();
            permColors[[1]] = permRefColors;
            permColors[[2]] = NA;
            permColorsForAcc = list();
            permColorsForAcc[[1]] = colorRef;
            permColorsForAcc[[2]] = colorTest;

            # For reproducibility of previous results, write separate code for threaded and unthreaded
            # calculations

            if (parallelCalculation)
            {
               combineCalculations = function(...)
               {
                 list(...);
               }
               seed = sample(1e8, 1);
               if (verbose > 2) 
                  printFlush(paste(spaces, " ......parallel calculation of permuted statistics.."));
               datout = foreach(perm = 1:nPermutations, .combine = combineCalculations,
                                .multicombine = TRUE, .maxcombine = nPermutations+10)%dopar% 
               {
                      set.seed(seed + perm + perm^2); 
                      collectGarbage();
                      .modulePreservationInternal(permExpr, permColors, dataIsExpr = dataIsExpr,
                                                  calculatePermutation = TRUE,
                                                  multiColorForAccuracy = permColorsForAcc,
                                                  networkType = networkType,
                                                  corFnc = corFnc, corOptions = corOptions,
                                                  referenceNetworks = 1, 
                                                  testNetworks = list(2),
                                                  densityOnly = useInterpolation,
                                                  maxGoldModuleSize = maxGoldModuleSize,
                                                  maxModuleSize = maxModuleSize, quickCor = quickCor,
                                                  ccTupletSize = ccTupletSize,
                                                  calculateCor.kIMall = calculateCor.kIMall,
                                                  calculateClusterCoeff = calculateClusterCoeff,
                                                  # calculateQuality = calculateQuality,
                                                  greyName = greyName,
                                                  checkData = FALSE,
                                                  verbose = verbose -3, indent = indent + 3)
               }
               for (perm in 1:nPermutations)
               {
                 if (!datout[[perm]] [[1]]$netPresent[2])
                 stop(paste("Internal error: no data in permuted set preservation measures. \n",
                            "Please contact the package maintainers. Sorry!"))
                 permOut[[iref]][[tnet]]$regStats[, , perm] = as.matrix(
                           cbind(datout[[perm]] [[1]]$quality[[2]][, -1],
                                 datout[[perm]] [[1]]$intra[[2]],
                                 datout[[perm]] [[1]]$inter[[2]]));
                 permOut[[iref]][[tnet]]$fixStats[, , perm] = as.matrix(datout[[perm]] [[1]]$accuracy[[2]]);
               }
               datout = datout[[1]] # For the name stting procedures that follow...
            } else for (perm in 1:nPermutations ) 
            {
               if (verbose > 2) printFlush(paste(spaces, " ......working on permutation", perm));
               #newRNG = .Random.seed;
               #if (!is.null(oldRNG))
               #  if (isTRUE(all.equal(newRNG, oldRNG)))
               #    printFlush("WARNING: something's wrong with the RNG... old and new RNG equal.");
               #oldRNG = .Random.seed
               #set.seed(perm*2);
               datout= .modulePreservationInternal(permExpr, permColors, dataIsExpr = dataIsExpr, 
                                                  calculatePermutation = TRUE,
                                                  multiColorForAccuracy = permColorsForAcc, 
                                                  networkType = networkType,
                                                  corFnc = corFnc, corOptions = corOptions,
                                                  referenceNetworks=1, 
                                                  testNetworks = list(2),
                                                  densityOnly = useInterpolation,
                                                  maxGoldModuleSize = maxGoldModuleSize,
                                                  maxModuleSize = maxModuleSize, quickCor = quickCor,
                                                  ccTupletSize = ccTupletSize,
                                                  calculateCor.kIMall = calculateCor.kIMall,
                                                  calculateClusterCoeff = calculateClusterCoeff,
                                                  # calculateQuality = calculateQuality,
                                                  greyName = greyName,
                                                  checkData = FALSE, 
                                                  verbose = verbose -3, indent = indent + 3)
               if (!datout[[1]]$netPresent[2])
                  stop(paste("Internal error: no data in permuted set preservation measures. \n",
                             "Please contact the package maintainers. Sorry!"))
               permOut[[iref]][[tnet]]$regStats[, , perm] = as.matrix(cbind(datout[[1]]$quality[[2]][, -1], 
                                                                  datout[[1]]$intra[[2]], 
                                                                  datout[[1]]$inter[[2]]));
               permOut[[iref]][[tnet]]$fixStats[, , perm] = as.matrix(datout[[1]]$accuracy[[2]]);
               collectGarbage();
            }
            regStatNames = c(colnames(datout[[1]]$quality[[2]])[-1], colnames(datout[[1]]$intra[[2]]), 
                             colnames(datout[[1]]$inter[[2]]));
            regModuleNames[[iref]][[tnet]] = rownames(datout[[1]]$quality[[2]]);
            regModuleSizes[[iref]][[tnet]] = datout[[1]]$quality[[2]][, 1]
            fixStatNames = colnames(datout[[1]]$accuracy[[2]]);
            fixModuleNames = rownames(datout[[1]]$accuracy[[2]]);
            dimnames(permOut[[iref]][[tnet]]$regStats) = list(regModuleNames[[iref]][[tnet]], regStatNames,
                                                             spaste("Permutation.", c(1:nPermutations)));
            dimnames(permOut[[iref]][[tnet]]$fixStats) = list(fixModuleNames, fixStatNames,
                                                             spaste("Permutation.", c(1:nPermutations)));
            permutationsPresent[tnet, iref] = TRUE
          } else {
            regModuleNames[[iref]][[tnet]] = NA;
            regModuleSizes[[iref]][[tnet]] = NA;
            permOut[[iref]][[tnet]] = NA;
          }
      }
         
      if (savePermutedStatistics) 
         save(regModuleSizes, regStatNames, regModuleNames, fixStatNames, 
              fixModuleNames, permutationsPresent, interpolationUsed, permOut, file=permutedStatisticsFile)
   }  # if (!psLoaded)

   collectGarbage();

   if (any(interpolationUsed, na.rm = TRUE))
   {
      if (verbose > 0) printFlush(paste(spaces, "..Calculating interpolation approximations.."));
      if (plotInterpolation) pdf(file = interpolationPlotFile, width = 12, height = 6)
   }

   # Define a "class" indicating that no valid fit was obtained

   invalidFit = "invalidFit";
     
   LogIndex = c(rep(c(1, 0, 0, 0, 0, 0, 0), 2), 0, 0, 0, 0, 0, 0)
   dfMean = c( rep(c(2, 2, 2, 1, 1, 1, 1), 2), 3, 3, 3, 2, 2, 2)
   dfSD = c( rep(c(2, 2, 2, 2, 2, 2, 2), 2), 2, 2, 2, 2, 2, 2)

   epsilon=0.00001
   meanLM=list()
   seLM=list()
   OmitModule=c(permGoldName,permGreyName)
   for (iref in 1:nRefNets) 
   {
      ref = referenceNetworks[iref];
      meanLM[[iref]]=list()
      seLM[[iref]]=list()
      for (tnet in 1:nNets) if (observed[[iref]]$netPresent[tnet] & 
                                permutationsPresent[tnet, iref] & interpolationUsed[tnet, iref]) 
      {
         NotGold = !is.element(regModuleNames[[iref]][[tnet]], OmitModule)
         meanLM[[iref]][[tnet]] = list()
         seLM[[iref]][[tnet]] = list()
         logName=matrix("", sum(NotGold), 1) 
         for (stat in 1:nRegStats)
         {
            means = c(apply(permOut[[iref]][[tnet]]$regStats[NotGold, stat, , drop = FALSE], c(1:2), 
                         mean, na.rm=TRUE));
            SD=try( c(apply(permOut[[iref]][[tnet]]$regStats[NotGold, stat, , drop = FALSE], c(1:2), 
                         sd, na.rm = TRUE)),
                    silent = TRUE);
            if (class(SD)=='try-error') SD = NA;
            if (any(is.na(c(means, SD))) )
            {
               meanLM[[iref]][[tnet]][[stat]] = NA;
               class(meanLM[[iref]][[tnet]][[stat]]) = invalidFit;
               seLM[[iref]][[tnet]][[stat]] = NA;
               class(seLM[[iref]][[tnet]][[stat]]) = invalidFit;
            } else {
               modSizes = regModuleSizes[[iref]][[tnet]][NotGold];
               xx = log(modSizes)
               if(LogIndex[stat]==1) 
               {
                  means[means<epsilon]=epsilon
                  yy=log(means)
                  logName[stat] = "log"
               } else{
                  yy=means
                  logName[stat] = ""
               }
               name1=regStatNames[stat];
               name2 = paste(setNames[ref], "vs", setNames[tnet]);
               meanLM[[iref]][[tnet]][[stat]]=lm(yy~ ns(xx, df=dfMean[stat] )) 
               names( meanLM[[iref]][[tnet]])[stat]=paste(name1,"_meanInterpolation",sep="")
               PredictedMedian=as.numeric( predict(meanLM[[iref]][[tnet]][[stat]]))
               if(all(!is.na(means))&&length(table(means))>1)
               {
                  par(mfrow=c(1,2))
                  par(mar = c(5, 5, 4, 1));
                  SE = SD/sqrt(nPermutations);
                  ymin = min(yy-SE, na.rm = TRUE)
                  ymax = max(yy+SE, na.rm = TRUE)
                  verboseScatterplot(xx,yy,xlab="Log module size",
                                     ylab=paste(logName[stat],"Permutation median"),
                                     main = paste("mean", name1, "\n", name2, "; "),
                                     cex.axis = 1, cex.lab = 1, cex.main = 1.2,
                                     ylim = c(ymin, ymax))
                  try(errbar(xx,yy,yy+SE, yy-SE, add = TRUE), silent = TRUE)
                  lines(xx[order(xx)], PredictedMedian[order(xx)] , col="red")
                  verboseScatterplot(yy, PredictedMedian,xlab="Observed permutation median",
                                     ylab="Predicted permutation median",
                                     main = paste("mean", name1, "\n", name2, "\n"),
                                     cex.axis = 1, cex.lab = 1, cex.main = 1.2)
                  abline(0,1,col="green")
               }
               epsilon2=0.0000000001
               SD[SD<epsilon2]=epsilon2
                  
               yy2=log(SD)
               xx2=log(modSizes)
               seLM[[iref]][[tnet]][[stat]]=lm(yy2~ ns(xx2, df=dfSD[stat] ) )
               names(seLM[[iref]][[tnet]])[stat]=paste(name1,"_SDInterpolation",sep="")
               PredictedSD=as.numeric( predict(seLM[[iref]][[tnet]][[stat]]))
               if(all(!is.na(SD))&&length(table(SD))>1)
               {
                     par(mfrow=c(1,2))
                     verboseScatterplot(xx2,yy2,xlab="Log Module Size",ylab="Log observed permutation SD",
                                     main = paste("SD", name1, "\n", name2, "\n"),
                                     cex.axis = 1, cex.lab = 1, cex.main = 1.2)
                     lines(xx2[order(xx2)], PredictedSD[order(xx2)] , col="red")
                     verboseScatterplot(yy2, PredictedSD,xlab="Log observed permutation SD",
                                        ylab="Log predicted permutation SD",
                                     main = paste("SD", name1, "\n", name2, "\n"),
                                     cex.axis = 1, cex.lab = 1, cex.main = 1.2)
                     abline(0,1,col="green")
               }
            } 
         }
      }      
   }
       
   if (any(interpolationUsed))
   {
      if (plotInterpolation) dev.off();
   }

   observedQuality = list();
   observedPreservation = list();
   observedReferenceSeparability = list();
   observedTestSeparability = list();
   observedOverlapCounts = list();
   observedOverlapPvalues = list();
   observedAccuracy = list();
   Z.quality = list();
   Z.preservation = list();
   Z.referenceSeparability = list();
   Z.testSeparability = list();
   Z.accuracy = list();
   interpolationStat = c(rep(TRUE, 14), FALSE, FALSE, FALSE, rep(TRUE, 6));
   for(iref in 1:nRefNets)
   {
     observedQuality[[iref]] = list();
     observedPreservation[[iref]] = list();
     observedReferenceSeparability[[iref]] = list();
     observedTestSeparability[[iref]] = list();
     observedOverlapCounts[[iref]] = list();
     observedOverlapPvalues[[iref]] = list();
     observedAccuracy[[iref]] = list();
     Z.quality[[iref]] = list();
     Z.preservation[[iref]] = list();
     Z.referenceSeparability[[iref]] = list();
     Z.testSeparability[[iref]] = list();
     Z.accuracy[[iref]] = list();
     for (tnet in 1:nNets) if (observed[[iref]]$netPresent[tnet] & permutationsPresent[tnet, iref])
     {
       nModules = nrow(observed[[iref]]$intra[[tnet]]);
       nQualiStats = ncol(observed[[iref]]$quality[[tnet]])-1;
       nIntraStats = ncol(observed[[iref]]$intra[[tnet]]);
       accuracy = observed[[iref]]$accuracy[[tnet]]
       inter = observed[[iref]]$inter[[tnet]];
       nInterStats = ncol(inter) + ncol(accuracy);
       a2i = match(rownames(accuracy), rownames(inter));
       accuracy2 = matrix(NA, nrow(inter), ncol(accuracy));
       accuracy2[a2i, ] = accuracy;
       rownames(accuracy2) = rownames(inter);
       colnames(accuracy2) = colnames(accuracy);
       modSizes = observed[[iref]]$quality[[tnet]][, 1]
       allObsStats = cbind(observed[[iref]]$quality[[tnet]][, -1], 
                           observed[[iref]]$intra[[tnet]], 
                           accuracy2, inter);
       sepCol = match("separability.qual", colnames(observed[[iref]]$quality[[tnet]]));
       sepCol2 = match("separability.pres", colnames(observed[[iref]]$intra[[tnet]]));
       quality = observed[[iref]]$quality[[tnet]][, -sepCol];
       if (dataIsExpr | (!restrictSummaryForGeneralNetworks))
       {
         rankColsQuality = c(2,3,4,5)
         rankColsDensity = c(1,2,3,4)
       } else {
         rankColsQuality = c(5)
         rankColsDensity = 4;
       }
     
       ranks = apply(-quality[, rankColsQuality, drop = FALSE], 2, rank, na.last = "keep");
       medRank = apply(as.matrix(ranks), 1, median, na.rm = TRUE);
       observedQuality[[iref]][[tnet]] = cbind(moduleSize = modSizes,
                                               medianRank.qual = medRank,
                                               quality[, -1]);
       preservation = cbind(observed[[iref]]$intra[[tnet]][, -sepCol2], inter );
    
       ranksDensity = apply(-preservation[, rankColsDensity, drop = FALSE], 2, rank, na.last = "keep");
       medRankDensity = apply(as.matrix(ranksDensity), 1, median, na.rm = TRUE);
       if (dataIsExpr | (!restrictSummaryForGeneralNetworks))
       {
         if (includekMEallInSummary)
         {
            connSummaryInd = c(7:10)
         } else {
            connSummaryInd = c(7,8,10);
         }
       } else {
         connSummaryInd = c(7,10) # in this case only cor.kIM and cor.Adj which sits in the cor.cor slot
       }
       ranksConnectivity = apply(-preservation[, connSummaryInd, drop = FALSE], 2, rank, na.last = "keep");
       medRankConnectivity = apply(as.matrix(ranksConnectivity), 1, median, na.rm = TRUE);
       medRank = apply(cbind(ranksDensity, ranksConnectivity), 1, median, na.rm = TRUE);
       observedPreservation[[iref]][[tnet]] = cbind(moduleSize = modSizes, 
                                                    medianRank.pres = medRank,
                                                    medianRankDensity.pres = medRankDensity,
                                                    medianRankConnectivity.pres = medRankConnectivity,
                                                    preservation);
       observedAccuracy[[iref]][[tnet]] = cbind(moduleSize = modSizes[a2i], accuracy);
       observedReferenceSeparability[[iref]][[tnet]] = cbind(moduleSize = modSizes, 
                                                     observed[[iref]]$quality[[tnet]][sepCol]);
       observedTestSeparability[[iref]][[tnet]] = cbind(moduleSize = modSizes, 
                                                observed[[iref]]$intra[[tnet]][sepCol2]);
       
       observedOverlapCounts[[iref]][[tnet]] = observed[[iref]]$overlapTables[[tnet]]$countTable;
       observedOverlapPvalues[[iref]][[tnet]] = observed[[iref]]$overlapTables[[tnet]]$pTable;

       nAllStats = ncol(allObsStats);
       zAll = matrix(NA, nModules, nAllStats);
       rownames(zAll) = rownames(observed[[iref]]$intra[[tnet]])
       colnames(zAll) = paste("Z.", colnames(allObsStats), sep="");
       logModSizes = log(modSizes);
       goldRowPerm = match(goldName, regModuleNames[[iref]][[tnet]]);
       goldRowObs = match(goldName, rownames(inter));
       fixInd = 1;
       regInd = 1;
       for (stat in 1:nAllStats) 
         if (interpolationStat[stat])
         {
           if (interpolationUsed[tnet, iref])
           {
              if (class(meanLM[[iref]][[tnet]][[regInd]])!=invalidFit)
              {
                 #print(paste(regInd, stat))
                 prediction = predict(meanLM[[iref]][[tnet]][[regInd]],
                                           newdata = data.frame(xx = logModSizes), se.fit = TRUE)
                 predictedMean = as.numeric(prediction$fit);
                 if (LogIndex[regInd]==1) 
                    predictedMean = exp(predictedMean);
                 predictedSD = exp(as.numeric(predict(seLM[[iref]][[tnet]][[regInd]],
                                                   newdata = data.frame(xx2 = logModSizes)))); 
                 zAll[, stat] = (allObsStats[, stat] - predictedMean)/predictedSD
                 # For the gold module : take the direct observations.
                 goldMean = mean(permOut[[iref]][[tnet]]$regStats[goldRowPerm, regInd, ], na.rm = TRUE);
                 goldSD = apply(permOut[[iref]][[tnet]]$regStats[goldRowPerm, regInd, ], 2, sd, na.rm = TRUE);
                 zAll[goldRowObs, stat] = (allObsStats[goldRowObs, stat] - goldMean)/goldSD
              }
           } else {
              means = c(apply(permOut[[iref]][[tnet]]$regStats[, regInd, , drop = FALSE], 
                        c(1:2), mean, na.rm = TRUE));
              SDs = c(apply(permOut[[iref]][[tnet]]$regStats[, regInd, , drop = FALSE], 
                        c(1:2), sd, na.rm = TRUE));
              z = ( allObsStats[, stat] - means) / SDs;
              if (any(is.finite(z)))
              {
                 finite = is.finite(z)
                 z[finite][SDs[finite]==0] = max(abs(z[finite]), na.rm = TRUE) *
                               sign(allObsStats[, stat] - means)[SDs[finite]==0]
              }
              zAll[, stat] = z;
           }
           regInd = regInd + 1;
         } else {
            if (.cvPresent(multiColor[[tnet]]) )
            {
               means = c(apply(permOut[[iref]][[tnet]]$fixStats[, fixInd, , drop = FALSE], 
                              c(1:2), mean, na.rm = TRUE));
               SDs = c(apply(permOut[[iref]][[tnet]]$fixStats[, fixInd, , drop = FALSE], 
                              c(1:2), sd, na.rm = TRUE));
               z = ( allObsStats[a2i, stat] - means) / SDs;
               if (any(is.finite(z)))
               {
                  finite = is.finite(z)
                  z[finite][SDs[finite]==0] = max(abs(z[finite]), na.rm = TRUE) *
                                sign(allObsStats[a2i, stat] - means)[SDs[finite]==0]
               }
               zAll[a2i, stat] = z;
            }
            fixInd = fixInd + 1
         }
       zAll = as.data.frame(zAll);
       sepCol = match("Z.separability.qual", colnames(zAll)[1:nQualiStats]);
       zQual = zAll[, c(1:nQualiStats)][ , -sepCol];
       summaryColsQuality = rankColsQuality -1 # quality also contains module sizes, Z does not
       summZ = apply(zQual[, summaryColsQuality, drop = FALSE], 1, median, na.rm = TRUE); 
       Z.quality[[iref]][[tnet]] = data.frame(cbind(moduleSize = modSizes, Zsummary.qual = summZ, 
                                                    zQual));
       Z.referenceSeparability[[iref]][[tnet]] = 
             data.frame(cbind(moduleSize = modSizes, zAll[sepCol]));
       st = nQualiStats + 1; en = nQualiStats + nIntraStats + nInterStats;
       sepCol2 = match("Z.separability.pres", colnames(zAll)[st:en]);
       accuracyCols = match(c("Z.accuracy", "Z.minusLogFisherP", "Z.coClustering"), colnames(zAll)[st:en]);
       zPres = zAll[, st:en][, -c(sepCol2, accuracyCols)];
       summaryColsPreservation = list(summaryColsQuality, connSummaryInd);
       nGroups = length(summaryColsPreservation)
       summZMat = matrix(0, nrow(zPres), nGroups);
       for (g in 1:nGroups)
         summZMat[, g] = apply(zPres[, summaryColsPreservation[[g]], drop = FALSE], 1, median, na.rm = TRUE);
       colnames(summZMat) = c("Zdensity.pres", "Zconnectivity.pres");
       summZ = apply(summZMat, 1, mean, na.rm = TRUE);
       Z.preservation[[iref]][[tnet]] = data.frame(cbind(moduleSize = modSizes, Zsummary.pres = summZ, 
                                                         summZMat, zPres));
       Z.testSeparability[[iref]][[tnet]] = 
             data.frame(cbind(moduleSize = modSizes, zAll[nQualiStats + sepCol2]));
       Z.accuracy[[iref]][[tnet]] = 
             data.frame(cbind(moduleSize = modSizes, zAll[st:en][accuracyCols]));
     } else {
       observedQuality[[iref]][[tnet]] = NA;
       observedPreservation[[iref]][[tnet]] = NA;
       observedReferenceSeparability[[iref]][[tnet]] = NA;
       observedTestSeparability[[iref]][[tnet]] = NA;
       observedAccuracy[[iref]][[tnet]] = NA;
       observedOverlapCounts[[iref]][[tnet]] = NA;
       observedOverlapPvalues[[iref]][[tnet]] = NA;
       Z.quality[[iref]][[tnet]] = NA;
       Z.preservation[[iref]][[tnet]]= NA;
       Z.referenceSeparability[[iref]][[tnet]]= NA;
       Z.testSeparability[[iref]][[tnet]]= NA;
       Z.accuracy[[iref]][[tnet]] = NA;
     }
     names(observedQuality[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedPreservation[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedReferenceSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedTestSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedAccuracy[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedOverlapCounts[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedOverlapPvalues[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(Z.quality[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                      setNames);
     names(Z.preservation[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                      setNames);
     names(Z.referenceSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                                    setNames);
     names(Z.testSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                                    setNames);
     names(Z.accuracy[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                                    setNames);
   } 

   names(observedQuality) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedPreservation) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedReferenceSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedTestSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedAccuracy) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedOverlapCounts) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedOverlapPvalues) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.quality) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.preservation) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.referenceSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.testSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.accuracy) = paste("ref", setNames[referenceNetworks], sep=".");

   summaryIndQuality = list(c(3:6));
   summaryIndPreservation = list(c(5:8), connSummaryInd + 4, c(3,4)); #4 = 1 module size + 3 summary indices
   p.quality = .pValueFromZ(Z.quality, summaryCols = 2, summaryInd = summaryIndQuality);
   p.preservation = .pValueFromZ(Z.preservation, summaryCols = c(3,4,2), summaryInd = summaryIndPreservation);
   p.referenceSeparability = .pValueFromZ(Z.referenceSeparability);
   p.testSeparability = .pValueFromZ(Z.testSeparability);
   p.accuracy = .pValueFromZ(Z.accuracy);

   pBonf.quality = .pValueFromZ(Z.quality, bonf = TRUE, summaryCols = 2, summaryInd = summaryIndQuality);
   pBonf.preservation = .pValueFromZ(Z.preservation, bonf = TRUE, summaryCols = c(3,4,2),
                                     summaryInd = summaryIndPreservation)
   pBonf.referenceSeparability = .pValueFromZ(Z.referenceSeparability, bonf = TRUE);
   pBonf.testSeparability = .pValueFromZ(Z.testSeparability, bonf = TRUE);
   pBonf.accuracy = .pValueFromZ(Z.accuracy, bonf = TRUE);

   if (calculateQvalue)
   {
     q.quality = .qValueFromP(p.quality, summaryCols = 2, summaryInd = summaryIndQuality)
     q.preservation = .qValueFromP(p.preservation, summaryCols = c(3,4,2), 
                                   summaryInd = summaryIndPreservation);
     q.referenceSeparability = .qValueFromP(p.referenceSeparability);
     q.testSeparability = .qValueFromP(p.testSeparability);
     q.accuracy = .qValueFromP(p.accuracy);
   } else {
     q.quality = NULL;
     q.preservation = NULL;
     q.referenceSeparability = NULL;
     q.testSeparability = NULL;
     q.accuracy = NULL;
   }

   output=list(quality = list(observed = observedQuality, Z = Z.quality, 
                              log.p = p.quality, 
                              log.pBonf = pBonf.quality,
                              q = q.quality ),
               preservation = list(observed = observedPreservation, Z = Z.preservation,
                                   log.p = p.preservation,
                                   log.pBonf = pBonf.preservation,
                                   q = q.preservation),
               accuracy = list(observed = observedAccuracy, Z = Z.accuracy,
                               log.p = p.accuracy,
                               log.pBonf = pBonf.accuracy,
                               q = q.accuracy,
                               observedCounts = observedOverlapCounts,
                               observedFisherPvalues = observedOverlapPvalues),
               referenceSeparability = list(observed = observedReferenceSeparability,
                                            Z = Z.referenceSeparability,
                                            log.p = p.referenceSeparability,
                                            log.pBonf = pBonf.referenceSeparability,
                                            q = q.referenceSeparability),
               testSeparability = list(observed = observedTestSeparability,
                                       Z = Z.testSeparability,
                                       log.p = p.testSeparability,
                                       log.pBonf = pBonf.testSeparability,
                                       q = q.testSeparability),
               permutationDetails = list(permutedStatistics = permOut, 
                                         interpolationModuleSizes = regModuleSizes,
                                         interpolationStatNames = regStatNames,
                                         permutationsPresent = permutationsPresent,
                                         interpolationUsed = interpolationUsed)
              );

   checkComps = list(c(1,2), c(1,2), c(1,2), c(1,2), c(1,2));
   nCheckComps = length(checkComps);
   if (discardInvalidOutput)
   {
     for(iref in 1:nRefNets)
     {
       for (tnet in 1:nNets) if (observed[[iref]]$netPresent[tnet] & permutationsPresent[tnet, iref])
       {
         for (oc in 1:nCheckComps)
         {
            for (ic in checkComps[[oc]])
            {
               data = output[[oc]][[ic]][[iref]][[tnet]]
               keep = apply(!is.na(data), 2, sum) > 0;
               output[[oc]][[ic]][[iref]][[tnet]] = data[, keep, drop = FALSE];
            }
         }
       }
     } 
   }

   return(output)
}


#=====================================================================================================
#
# .modulePreservationInternal
#
#=====================================================================================================

# Calculate module preservation scores for a given multi-expression data set.


# color vector present?

.cvPresent = function(cv)
{
  if (is.null(cv)) return(FALSE);
  if (length(cv)==1 && (is.na(cv[1]))) return(FALSE);
  return(TRUE);
}

.nNAColors = function(cv) {if (.cvPresent(cv)) sum(is.na(cv)) else 0}


.accuracyStatistics = function(colorRef, colorTest, ccTupletSize, greyName, pEpsilon)
{
   colorRefLevels = sort(unique(colorRef));
   nRefMods = length(colorRefLevels);
   refModSizes = table(colorRef);

   accuracy=matrix(NA,nRefMods ,3)
   colnames(accuracy)=c("accuracy", "minusLogFisherP", "coClustering");
   rownames(accuracy)=colorRefLevels;

   nRefGenes = length(colorRef); # also equals nTestGenes
   if(.cvPresent(colorTest))
   {
      #if (verbose > 1) printFlush(paste(spaces, "....calculating color label accuracy..."));
      overlap = overlapTable(colorTest, colorRef)
      greyRow = rownames(overlap$countTable)==greyName;
      greyCol = colnames(overlap$countTable)==greyName;
      #bestCount = apply(overlap$countTable, 2, max);
      overlap$pTable[!is.finite(overlap$pTable)] = pEpsilon;
      overlap$pTable[overlap$pTable < pEpsilon] = pEpsilon;
      bestCount = rep(0, ncol(overlap$countTable));
      if (sum(!greyRow) > 0 & sum(!greyCol)>0)
      {
         bestCount[!greyCol] = apply(overlap$countTable[!greyRow, !greyCol, drop = FALSE], 2, max)
         accuracy[!greyCol, 2] = 
                           -log10(apply(overlap$pTable[!greyRow, !greyCol, drop = FALSE], 2, min));
         ccNumer = apply(overlap$countTable[!greyRow, !greyCol, drop = FALSE], 2, choose,
                         ccTupletSize);
         ccDenom = choose(refModSizes[!greyCol], ccTupletSize)
         accuracy[!greyCol, 3] = apply(ccNumer, 2, sum)/ccDenom;
      }
      if (sum(greyRow)==1 & sum(greyCol) ==1)
      {
         bestCount[greyCol] = overlap$countTable[greyRow, greyCol];
         accuracy[greyCol, 2] = -log10(overlap$pTable[greyRow, greyCol]);
         accuracy[greyCol, 3] = choose(bestCount[greyCol], ccTupletSize)/
                                  choose(refModSizes[greyCol], ccTupletSize);
      }
      accuracy[, 1] = bestCount/refModSizes;
   } else {
      overlap = list(countTable = NA, pTable = NA);
   }
   list(accuracy = accuracy, overlapTable = overlap);
}

#=================================================================================================
#
# .modulePreservationInternal
#
#=================================================================================================

# multiData contains either expression data or adjacencies; which one is indicated in dataIsExpr

.modulePreservationInternal = function(multiData, multiColor, dataIsExpr,
                              calculatePermutation,
                              multiColorForAccuracy = NULL,
                              networkType = "signed", 
                              corFnc = "cor", 
                              corOptions = "use = 'p'",
                              referenceNetworks=1,
                              testNetworks,
                              densityOnly = FALSE,
                              maxGoldModuleSize = 1000, 
                              maxModuleSize = 1000, quickCor = 1,
                              ccTupletSize,
                              calculateCor.kIMall,
                              calculateClusterCoeff,
                              # calculateQuality = FALSE, 
                              checkData = TRUE,
                              greyName,
                              pEpsilon = 1e-200,
                              verbose = 1, indent = 0)
{

   spaces = indentSpaces(indent);

   #size = checkSets(multiData);

   nNets = length(multiData);
   nGenes = sapply(multiData, sapply, ncol);

   nType = charmatch(networkType, .networkTypes);
   if (is.na(networkType))
      stop(paste("Unrecognized networkType argument.",
           "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")));

   setNames = names(multiData);
   # Check multiColor. 

   if (length(multiColor)!=nNets)
   {
     multiColor2=list()
     if (length(names(multiColor))!=length(multiColor))
       stop("Each entry of 'multiColor' must have a name.");
     color2expr = match(names(multiColor), names(multiData));
     if (any(is.na(color2expr)))
       stop("Entries of 'multiColor' must name-match entries in 'multiData'.");
     for(s in 1:nNets)
     {
        multiData[[s]]$data = as.matrix(multiData[[s]]$data);
        loc = which(names(multiColor) %in% names(multiData)[s])
        if (length(loc)==0) 
        {
          multiColor2[[s]] = NA
        } else 
          multiColor2[[s]]=multiColor[[loc]]
     }
     multiColor = multiColor2
     rm(multiColor2);
   }

   if (is.numeric(greyName)) goldName = 0.1 else goldName = "gold"

   MEgrey = paste("ME", greyName, sep="");
   MEgold = paste("ME", goldName, sep="");
     
   collectGarbage();

   for (s in 1:nNets)
   {
     if ( .cvPresent(multiColor[[s]]) & nGenes[s]!=length(multiColor[[s]]))
       stop(paste("Color vector for set", s, "does not have the correct number of entries."));
   }

   keepGenes = list();
   for (s in 1:nNets) 
      keepGenes[[s]] = rep(TRUE, nGenes[s])

   nNAs = sapply(multiColor, .nNAColors)
   if (any(nNAs > 0))
   {
     if (verbose > 0) printFlush(paste(spaces, " ..removing genes with missing color..."))
     for (s in 1:nNets) 
       if (.cvPresent(multiColor[[s]]))
       {
          keepGenes[[s]] = !is.na(multiColor[[s]]);
          if (dataIsExpr)
          {
             multiData[[s]]$data = multiData[[s]]$data[, keepGenes[[s]]];
          } else {
             multiData[[s]]$data = multiData[[s]]$data[keepGenes[[s]], keepGenes[[s]]];
          }
          multiColor[[s]] = multiColor[[s]][keepGenes[[s]]];
       } 
   }

   #if (verbose > 0) printFlush(paste(spaces, " ..preservation tests based on different reference networks"))
        
   datout=list()
   if (nNets==1) # && length(multiColor)==1)
   {     
      # For now stop; in the future we'll bring this up to speed as well.
      stop("Calculation of quality in an idividual network is not supported at this time. Sorry!");
   } else {  		
     for (iref in 1:length(referenceNetworks))
     {
       ref = referenceNetworks[iref]
       if (!.cvPresent(multiColor[[ref]]))
          stop(paste("Network", ref, 
                     "does not have a color vector and cannot be used as reference network."))
       if (verbose > 0) printFlush(paste(spaces, "..working on reference network",setNames[ref]))
       accuracy=list()
       quality = list();
       interPres=list()
       intraPres=list()
       overlapTables = list();

       netPresent = rep(FALSE, nNets)
       for (tnet in testNetworks[[iref]])
       {
          if (verbose > 1) printFlush(paste(spaces, "  ..working on test network",setNames[tnet]))
          overlap=intersect(colnames(multiData[[ref]]$data),colnames(multiData[[tnet]]$data))
          if (length(overlap)==0)
          {
            printFlush(paste(spaces, "WARNING: sets", ref, "and", tnet, 
                             "have no overlapping genes with valid colors.\n",
                             spaces, "No preservation measures can be calculated."));
            next;
          }
          loc1=match(overlap, colnames(multiData[[ref]]$data))
          loc2=match(overlap, colnames(multiData[[tnet]]$data))
          if(length(multiColor[[tnet]])>1)
          {
             colorTest=multiColor[[tnet]][loc2]
          } else  {
             colorTest=NA 
          }
          if (dataIsExpr)
          {
            datTest=multiData[[tnet]]$data[,loc2]
            datRef=multiData[[ref]]$data[,loc1]
          } else {
            datTest=multiData[[tnet]]$data[loc2, loc2]
            datRef=multiData[[ref]]$data[loc1, loc1]
          }
          colorRef=multiColor[[ref]][loc1]

          # if multiColorForAccuracy is present, use the colors in this list to calculate accuracy and
          # Fisher p values.

          if (!is.null(multiColorForAccuracy))
          {
            colorRefAcc = multiColorForAccuracy[[ref]][loc1];
            if (.cvPresent(multiColorForAccuracy[[tnet]]))
            {
               colorTestAcc = multiColorForAccuracy[[tnet]][loc2];
            } else 
               colorTestAcc = NA;
          } else {
            colorRefAcc = colorRef;
            colorTestAcc = colorTest; #This is the only place where colorTest is used (?)
          }

          # Accuracy measures

          if (calculatePermutation)
          {
             colorRefAcc = sample(colorRefAcc);
             colorTestAcc = sample(colorRefAcc);
          }
          x = .accuracyStatistics(colorRefAcc, colorTestAcc, 
                                 ccTupletSize = ccTupletSize, greyName = greyName, pEpsilon = pEpsilon);
          accuracy[[tnet]] = x$accuracy;
          overlapTables[[tnet]] = x$overlapTable;

          # From now on we work with colorRef; colorTest is not needed anymore.

          # Restrict each module to at most maxModuleSize genes..

          colorRefLevels = sort(unique(colorRef));
          nRefMods = length(colorRefLevels);
          nRefGenes = length(colorRef);

          # Check that the gold module is not too big. In particular, the gold module must not contain all valid
          # genes, since in such a case the random sampling makes no sense for density-based statistics.

          goldModSize = maxGoldModuleSize
          if (goldModSize > nRefGenes/2) goldModSize = nRefGenes/2;

          # ..step 1: gold module. Note that because of the above the gold module size is always smaller than
          # nRefGenes.

          goldModR = sample(nRefGenes, goldModSize)
          if (calculatePermutation)
          {
             goldModT = sample(nRefGenes, goldModSize)
          } else
             goldModT = goldModR
 
          if (dataIsExpr)
          {
            goldRef = datRef[, goldModR];
            goldRefP = datRef[, goldModT];
            goldTest = datTest[, goldModT];
          } else {
            goldRef = datRef[goldModR, goldModR];
            goldRefP = datRef[goldModT, goldModT];
            goldTest = datTest[goldModT, goldModT];
          }

          # ..step 2: proper modules and grey

          keepGenes = rep(TRUE, nRefGenes);
          for (m in 1:nRefMods)
          {
            inModule = colorRef == colorRefLevels[m]
            nInMod = sum(inModule)
            if(nInMod > maxModuleSize)
            {
               sam = sample(nInMod, maxModuleSize)
               keepGenes[inModule] = FALSE;
               keepGenes[inModule][sam] = TRUE;
            }
          }

          # Create the permuted data sets
          if (sum(keepGenes) < nRefGenes)
          {
            colorRef = colorRef[keepGenes]
             if (dataIsExpr)
             {
                datRef = datRef[, keepGenes]
             } else
                datRef = datRef[keepGenes, keepGenes]
             nRefGenes = length(colorRef);
             if (calculatePermutation)
             {
                keepPerm = sample(nRefGenes, sum(keepGenes));
                if (dataIsExpr)
                {
                   datTest = datTest[, keepPerm];
                   datRefP = datRef[, keepPerm];
                } else {
                   datTest = datTest[keepPerm, keepPerm];
                   datRefP = datRef[keepPerm, keepPerm];
                }
             } else {
                if (dataIsExpr)
                {
                  datTest = datTest[, keepGenes];
                  datRefP = datRef;
                } else {
                  datTest = datTest[keepGenes, keepGenes];
                  datRefP = datRef;
                }
             }
          } else {
             if (calculatePermutation)
             {
               perm = sample(c(1:nRefGenes));
               if (dataIsExpr)
               {
                 datRefP = datRef[, perm];
                 datTest = datTest[, perm];
               } else {
                 datRefP = datRef[perm, perm];
                 datTest = datTest[perm, perm];
               }
             } else {
               datRefP = datRef;
             }
          }
          if (dataIsExpr)
          {
             datRef = cbind(datRef, goldRef)
             datRefP = cbind(datRefP, goldRefP)
             datTest = cbind(datTest, goldTest)
          } else {
             datRef = .combineAdj(datRef, goldRef);
             datRefP = .combineAdj(datRefP, goldRefP);
             datTest = .combineAdj(datTest, goldTest);
             collectGarbage();
          }
             
          gold = rep(goldName, goldModSize)
          colorRef_2 = c(as.character(colorRef),gold)
          colorLevels = sort(unique(colorRef_2));
          opt = list(corFnc = corFnc, corOptions = corOptions, quickCor = quickCor, 
                     nType = nType, 
                     MEgold = MEgold, MEgrey = MEgrey, 
                     densityOnly = densityOnly, calculatePermutation = calculatePermutation,
                     calculateCor.kIMall = calculateCor.kIMall, 
                     calculateClusterCoeff = calculateClusterCoeff);
             
          if (dataIsExpr)
          {
            stats = .coreCalcForExpr(datRef, datRefP, datTest, colorRef_2, opt);
            interPresNames = spaste(corFnc, c(".kIM", ".kME", ".kMEall", 
                                             spaste(".", corFnc), ".clusterCoeff", ".MAR"));
            measureNames = c("propVarExplained", "meanSignAwareKME", "separability", 
                             "meanSignAwareCorDat", "meanAdj", "meanClusterCoeff", "meanMAR");
 
          } else {
            stats = .coreCalcForAdj(datRef, datRefP, datTest, colorRef_2, opt);
            interPresNames = spaste(corFnc, c(".kIM", ".kME", ".kIMall", ".adj", ".clusterCoeff", ".MAR"));
            measureNames = c("propVarExplained", "meanKIM", "separability", 
                             "meanSignAwareCorDat", "meanAdj", "meanClusterCoeff", "meanMAR");
          }

          name1=paste(setNames[[ref]],"_vs_",setNames[[tnet]],sep="")  
          quality[[tnet]] = cbind(stats$modSizes, 
                                  stats$proVar[, 1], 
                                  if (dataIsExpr) stats$meanSignAwareKME[, 1] else stats$meankIM[, 1],
                                  stats$Separability[, 1], stats$MeanSignAwareCorDat[,1],
                                  stats$MeanAdj[, 1], stats$meanClusterCoeff[, 1], 
                                  stats$meanMAR[, 1]);
          intraPres[[tnet]]=cbind(stats$proVar[, 2], 
                                  if (dataIsExpr) stats$meanSignAwareKME[, 2] else stats$meankIM[, 2],
                                  stats$Separability[, 2], stats$MeanSignAwareCorDat[, 2], 
                                  stats$MeanAdj[, 2], stats$meanClusterCoeff[, 2],
                                  stats$meanMAR[, 2])
          #colnames(quality[[tnet]]) = paste(c("moduleSize", measureNames), setNames[ref], sep = "_");
          #colnames(intraPres[[tnet]]) = paste(measureNames, name1, sep = "_");
          colnames(quality[[tnet]]) = c("moduleSize", paste(measureNames, "qual", sep="."));
          rownames(quality[[tnet]]) = colorLevels
          colnames(intraPres[[tnet]]) = paste(measureNames, "pres", sep=".");
          rownames(intraPres[[tnet]]) = colorLevels
          names(intraPres)[tnet]=paste(name1,sep="")
          quality[[tnet]] = as.data.frame(quality[[tnet]]);
          intraPres[[tnet]] = as.data.frame(intraPres[[tnet]]);
          interPres[[tnet]]= as.data.frame(cbind(stats$corkIM, stats$corkME, stats$corkMEall, stats$ICORdat,
                                                 stats$corCC, stats$corMAR))
          colnames(interPres[[tnet]])=interPresNames;
          rownames(interPres[[tnet]])=colorLevels
          names(interPres)[[tnet]]=paste(name1,sep="")
          netPresent[tnet] = TRUE;
      } # of for (test in testNetworks[[iref]])
              
      datout[[iref]]=list(netPresent = netPresent, quality = quality, 
                          intra = intraPres, inter = interPres, accuracy = accuracy,
                          overlapTables = overlapTables)
    } # of for (iref in 1:length(referenceNetworkss))
    names(datout)=setNames[referenceNetworks]
  }  # of else for if (nNets==1)
  return(datout)
        
}

.checkExpr = function(multiExpr, verbose, indent)
{
   spaces = indentSpaces(indent);
   nNets = length(multiExpr);
   if (verbose > 0) 
          printFlush(paste(spaces, " ..checking data for excessive amounts of missing data.."));
   for (set in 1:nNets)
   {
      gsg = goodSamplesGenes(multiExpr[[set]]$data, verbose= verbose -2, indent = indent + 2);
      if (!gsg$allOK)
      {
        stop(paste("The submitted 'multiExpr' data contain genes or samples\n",
              "  with zero variance or excessive counts of missing entries.\n",
              "  Please use the function goodSamplesGenes on each set to identify the problematic\n",
              "  genes and samples, and remove them before running modulePreservation."))
      }
   }
}

.checkAdj = function(multiAdj, verbose, indent)
{
   spaces = indentSpaces(indent);
   nNets = length(multiAdj);
   if (verbose > 0) 
     printFlush(paste(spaces, " ..checking adjacencies for excessive amounts of missing data"));
   for (set in 1:nNets)
   {
      checkAdjMat(multiAdj[[set]]$data);
      gsg = goodSamplesGenes(multiAdj[[set]]$data, verbose= verbose -2, indent = indent + 2);
      if (!gsg$allOK)
      {
        stop(paste("The submitted 'multiAdj' contains rows or columns\n",
              "  with zero variance or excessive counts of missing entries. Please remove\n",
              "  offending rows and columns before running modulePreservation."))
      }
   }
}

.combineAdj = function(block1, block2)
{
  n1 = ncol(block1);
  n2 = ncol(block2);
  comb = matrix(0, n1+n2, n1+n2);
  comb[1:n1, 1:n1] = block1;
  comb[(n1+1):(n1+n2), (n1+1):(n1+n2)] = block2;
  try( {colnames(comb) = c(colnames(block1), colnames(block2)) }, silent = TRUE);
  comb;
}


# This function is basically copied from the file networkConcepts.R

.computeLinksInNeighbors = function(x, imatrix){x %*% imatrix %*% x}
.computeSqDiagSum = function(x, vec) { sum(x^2 * vec) };

.clusterCoeff = function(adjmat1)
{
  # diag(adjmat1)=0
  no.nodes=dim(adjmat1)[[1]]
  nolinksNeighbors <- c(rep(-666,no.nodes))
  total.edge <- c(rep(-666,no.nodes))
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
  nolinksNeighbors <- apply(adjmat1, 1, .computeLinksInNeighbors, imatrix=adjmat1)
  subTerm = apply(adjmat1, 1, .computeSqDiagSum, vec = diag(adjmat1));
  plainsum  <- colSums(adjmat1)
  squaresum <- colSums(adjmat1^2)
  total.edge = plainsum^2 - squaresum
  #CChelp=rep(-666, no.nodes)
  CChelp=ifelse(total.edge==0,0, (nolinksNeighbors-subTerm)/total.edge)
  CChelp
} 

# This function assumes that the diagonal of the adjacency matrix is 1
.MAR = function(adjacency)
{
  denom = apply(adjacency, 2, sum)-1;
  mar = (apply(adjacency^2, 2, sum) - 1)/denom;
  mar[denom==0] = NA;
  mar;
}

    

#=======================================================================================
#
# Core calculation for expression data
#
#=======================================================================================

.coreCalcForExpr = function(datRef, datRefP, datTest, colors, opt)
{
  colorLevels = levels(factor(colors))
  nMods =length(colorLevels)
  nGenes = length(colors)

  # Flag modules whose size is 1
  modSizes = table(colors)
  act = (modSizes>1)

  kME=list()
          
  ME=list()
  if (!opt$densityOnly)
  {            
     ME[[1]]=moduleEigengenes(datRef, colors)$eigengenes
     kME[[1]]=as.matrix(signedKME(datRef,ME[[1]], corFnc = opt$corFnc, corOptions = opt$corOptions))
  }
   
  if (opt$calculatePermutation | opt$densityOnly)
  {
    #printFlush("Calculating ME[[2]]");
    ME[[2]]=moduleEigengenes(datRefP, colors)$eigengenes
    kME[[2]]=as.matrix(signedKME(datRefP,ME[[2]], corFnc = opt$corFnc, corOptions = opt$corOptions))
  } else {
    ME[[2]] = ME[[1]]
    kME[[2]] = kME[[1]]
  }

  ME[[3]]=moduleEigengenes(datTest, colors)$eigengenes
  kME[[3]]=as.matrix(signedKME(datTest,ME[[3]], corFnc = opt$corFnc, corOptions = opt$corOptions))

  modGenes = list();
  for (m in 1:nMods)
    modGenes[[m]] = c(1:nGenes)[colors==colorLevels[m]];
              
  #kME correlation
  #if (verbose > 1) printFlush(paste(spaces, "....calculating kME..."));
  corkME = rep(NA, nMods);
  corkMEall = rep(NA, nMods);
  if (!opt$densityOnly)
  {
     for(j in 1:nMods ) if(act[j])
     {
        nModGenes=length(modGenes[[j]]);
        names1=substring(colnames(kME[[1]]),4)
        j1=which(names1==colorLevels[j])
        #corkME[j]=cor(kME[[1]][loc,j1],kME[[3]][loc,j1], use = "p")
        corExpr = parse(text=paste(opt$corFnc, "(kME[[1]][modGenes[[j]],j1],kME[[3]][modGenes[[j]],j1]", 
                                   prepComma(opt$corOptions), ")"));
        corkME[j] = abs(eval(corExpr));

        #corkMEall[j]=cor(kME[[1]][,j1],kME[[3]][,j1], use = "p")
        corExpr = parse(text=paste(opt$corFnc, "(kME[[1]][,j1],kME[[3]][,j1] ", prepComma(opt$corOptions), ")"));
        corkMEall[j] = abs(eval(corExpr));
   #     covkME[j]=cov(kME[[1]][loc,j1],kME[[3]][loc,j1], use = "p")
   #     meanProductkME[j] = scalarProduct(kME[[1]][loc,j1],kME[[3]][loc,j1])
     }
  }
      
  #proportion of variance explained 
      
  #if (verbose > 1) 
  #  printFlush(paste(spaces, "....calculating proprotion of variance explained..."));
      
  proVar=matrix(NA, nMods ,2)
  meanSignAwareKME=matrix(NA, nMods ,2)
  names1=substring(colnames(kME[[2]]),4)
  for(j in 1:nMods ) if(act[j])
  {       
     j1=which(names1==colorLevels[j])
     proVar[j,1]=mean((kME[[2]][modGenes[[j]],j1])^2,na.rm=TRUE)
     proVar[j,2]=mean((kME[[3]][modGenes[[j]],j1])^2,na.rm=TRUE)
     if (opt$densityOnly)
     {
        if (opt$nType==1)
        {
          meanSignAwareKME[j,1]=mean(abs(kME[[2]][modGenes[[j]],j1]),na.rm = TRUE)
          meanSignAwareKME[j,2]=mean(abs(kME[[3]][modGenes[[j]],j1]),na.rm = TRUE)
        } else {
          meanSignAwareKME[j,1]=abs(mean(kME[[2]][modGenes[[j]],j1],na.rm = TRUE))
          meanSignAwareKME[j,2]=abs(mean(kME[[3]][modGenes[[j]],j1],na.rm = TRUE))
        }
     } else {
        meanSignAwareKME[j,1]=abs(mean(abs(kME[[2]][modGenes[[j]],j1]),na.rm = TRUE));
        meanSignAwareKME[j,2]=abs(mean(sign(kME[[1]][modGenes[[j]],j1]) * kME[[3]][modGenes[[j]],j1],na.rm =
TRUE))
     }
  }
  
  # if (verbose > 1) printFlush(paste(spaces, "....calculating separability..."));
  Separability=matrix(NA, nMods ,2)
  # for(k in (2-calculateQuality):2)
  for(k in 1:2)
  {
     Gold=which(colnames(ME[[k+1]]) %in% c(opt$MEgold, opt$MEgrey))
     #corME=cor(ME[[k+1]],use="p")
     corExpr = parse(text=paste(opt$corFnc, "(ME[[k+1]]", prepComma(opt$corOptions), ")"));
     corME= eval(corExpr);
     if (opt$nType==0) corME = abs(corME);
     diag(corME) = 0;
     Separability[,k]=1-apply(corME[, -Gold, drop = FALSE], 1, max, na.rm = TRUE)                        
  }

  #mean signed correlation&inter array correlation
  # if (verbose > 1) printFlush(paste(spaces, "....calculating MeanSignAwareCorDat..."));

  MeanSignAwareCorDat=matrix(NA,nMods ,2)
  ICORdat=rep(NA,nMods)
  corkIM = rep(NA, nMods);
  corCC = rep(NA, nMods);
  corMAR = rep(NA, nMods);
  MeanAdj = matrix(NA,nMods ,2)
  meanCC = matrix(NA,nMods ,2)
  meanMAR = matrix(NA,nMods ,2)
  for(j in 1:nMods ) if(act[j])
  {
     if (!opt$densityOnly) 
     {
        #ModuleCorData1=cor(datRef[,modGenes[[j]]],use="p", quick = as.numeric(opt$quickCor))
        corExpr = parse(text=paste(opt$corFnc, "(datRef[,modGenes[[j]]]", prepComma(opt$corOptions), 
                                           ", quick = as.numeric(opt$quickCor))"));
        ModuleCorData1=eval(corExpr);
     }
     if (opt$calculatePermutation | opt$densityOnly)
     {
        #ModuleCorData2=cor(datRefP[,modGenes[[j]]],use="p", quick = as.numeric(opt$quickCor))
        corExpr = parse(text=paste(opt$corFnc, "(datRefP[,modGenes[[j]]]", prepComma(opt$corOptions), 
                                           ", quick = as.numeric(opt$quickCor))"));
        ModuleCorData2 = eval(corExpr);
     } else 
        ModuleCorData2 = ModuleCorData1;

     #ModuleCorData3=cor(datTest[,modGenes[[j]]],use="p", quick = as.numeric(opt$quickCor))
     corExpr = parse(text=paste(opt$corFnc, "(datTest[,modGenes[[j]]]", prepComma(opt$corOptions),
                                           ", quick = as.numeric(opt$quickCor))"));
     ModuleCorData3 = eval(corExpr);
     if (opt$nType==1)
     {
        SignedModuleCorData2 = abs(ModuleCorData2)
     } else 
        SignedModuleCorData2 = ModuleCorData2;
     if (opt$densityOnly)
     {
        if (opt$nType==1) SignedModuleCorData3 = abs(ModuleCorData3)
             else SignedModuleCorData3 = ModuleCorData3
     } else
        SignedModuleCorData3 = sign(ModuleCorData1)*ModuleCorData3
     MeanSignAwareCorDat[j,1]=mean(as.dist(SignedModuleCorData2),na.rm = TRUE)
     MeanSignAwareCorDat[j,2]=mean(as.dist(SignedModuleCorData3),na.rm = TRUE)
     if (!opt$densityOnly)
     {
        #ICORdat[j]=cor(c(as.dist(ModuleCorData1)),c(as.dist(ModuleCorData3)),use="p")
        corExpr = parse(text=paste(opt$corFnc,
                                   "(c(as.dist(ModuleCorData1)),c(as.dist(ModuleCorData3))",
                                   prepComma(opt$corOptions), ")"));
        ICORdat[j] = eval(corExpr);
  #      ICOVdat[j]=cov(c(as.dist(ModuleCorData1)),c(as.dist(ModuleCorData3)),use="p")
  #      spdat[j]=scalarProduct(c(as.dist(ModuleCorData1)),c(as.dist(ModuleCorData3)))
     }

     if (opt$nType==1)
     {                            
        if (!opt$densityOnly) adjacency1 = ModuleCorData1^6;
        if (opt$calculatePermutation | opt$densityOnly) 
        {
            adjacency2 = ModuleCorData2^6;
        } else 
            adjacency2 = adjacency1;
        adjacency3 = ModuleCorData3^6;
     } else if (opt$nType==2)
     {
        if (!opt$densityOnly) adjacency1 = ( (1+ModuleCorData1)/2 ) ^12;
        if (opt$calculatePermutation | opt$densityOnly) 
        {
            adjacency2 = ( (1+ModuleCorData2)/2 ) ^12;
        } else 
            adjacency2 = adjacency1;
        adjacency3 = ( (1+ModuleCorData3)/2 ) ^12;
     } else {
        if (!opt$densityOnly) adjacency1 = ModuleCorData1^6;
        adjacency1[ModuleCorData1 < 0] = 0;
        if (opt$calculatePermutation | opt$densityOnly)
        {
           adjacency2 = ModuleCorData2^6;
           adjacency2[ModuleCorData2 < 0] = 0;
        } else
           adjacency2 = adjacency1;
        adjacency3 = ModuleCorData3^6;
        adjacency3[ModuleCorData3 < 0] = 0;
     }
     if (opt$calculateClusterCoeff)
     {
       ccRef = .clusterCoeff(adjacency1);
       ccRefP = .clusterCoeff(adjacency2);
       ccTest = .clusterCoeff(adjacency3);
       meanCC[j, 1] = mean(ccRefP);
       meanCC[j, 2] = mean(ccTest);
       if (!opt$densityOnly)
       {
         corExpr = parse(text=paste(opt$corFnc, "(ccRef, ccTest ", prepComma(opt$corOptions), ")"));
         corCC[j] = eval(corExpr);
       }
     } 
     marRef = .MAR(adjacency1);
     marRefP = .MAR(adjacency2);
     marTest = .MAR(adjacency3);
     meanMAR[j, 1] = mean(marRefP);
     meanMAR[j, 2] = mean(marTest);
     if (!opt$densityOnly)
     {
       kIMref = apply(adjacency1, 2, sum, na.rm = TRUE)
       kIMtest = apply(adjacency3, 2, sum, na.rm = TRUE)
       corExpr = parse(text=paste(opt$corFnc, "(kIMref, kIMtest ", prepComma(opt$corOptions), ")"));
       corkIM[j] = eval(corExpr);
       corExpr = parse(text=paste(opt$corFnc, "(marRef, marTest ", prepComma(opt$corOptions), ")"));
       corMAR[j] = eval(corExpr);
     }

     MeanAdj[j,1]=mean(as.dist(adjacency2), na.rm=TRUE)
     MeanAdj[j,2]=mean(as.dist(adjacency3), na.rm=TRUE)
  }  
  list(modSizes = modSizes, 
       corkIM = corkIM, corkME = corkME, corkMEall = corkMEall, 
       proVar = proVar, meanSignAwareKME = meanSignAwareKME,
       Separability = Separability, MeanSignAwareCorDat = MeanSignAwareCorDat, ICORdat = ICORdat,
       MeanAdj = MeanAdj, meanClusterCoeff = meanCC, meanMAR = meanMAR, corCC = corCC, corMAR = corMAR)
}

#===================================================================================================
#
# Core calculation for adjacency
#
#===================================================================================================

# A few supporting functions first:

.getSVDs = function(data, colors)
{
  colorLevels = levels(factor(colors))
  nMods =length(colorLevels)
  svds = list();
  for (m in 1:nMods)
  {
    modGenes = (colors==colorLevels[m])
    modAdj = data[modGenes, modGenes];
    if (sum(is.na(modAdj))>0)
    {
      seed = .Random.seed;
      modAdj = impute.knn(modAdj)$data;
      .Random.seed <<- seed;
    }
    svds[[m]] = svd(modAdj, nu=1, nv=0);
    svds[[m]]$u = c(svds[[m]]$u);
    if (sum(svds[[m]]$u, na.rm = TRUE) < 0) svds[[m]]$u = -svds[[m]]$u;
  }
  svds;
}

.kIM = function(adj, colors, calculateAll = TRUE)
{
  colorLevels = levels(factor(colors))
  nMods =length(colorLevels)
  nGenes = length(colors);
  kIM = matrix(NA, nGenes, nMods);
  if (calculateAll)
  {
     for (m in 1:nMods)
     {
       modGenes = colors==colorLevels[m];
       kIM[, m] = apply(adj[, modGenes, drop = FALSE], 1, sum, na.rm = TRUE);
       kIM[modGenes, m] = kIM[modGenes, m] - 1;
     }
  } else {
     for (m in 1:nMods)
     {
       modGenes = colors==colorLevels[m];
       kIM[modGenes, m] = apply(adj[modGenes, modGenes, drop = FALSE], 1, sum, na.rm = TRUE) - 1;
     }
  }
  kIM;
}


# Here is the main function

# Summary: 
#     PVE: from svd$d
#     kME: from svd$u (to make it different from kIM)
#     kIM: as usual
#     kMEall: from kIMall
#     meanSignAwareKME: from svd$u 
#     Separability: as in the paper
#     MeanSignAwareCorDat: ??
#     meanAdj: mean adjacency
#     cor.cor: replace by cor.adj

.coreCalcForAdj = function(datRef, datRefP, datTest, colors, opt)
{
#  printFlush(".coreCalcForAdj:entering");
  colorLevels = levels(factor(colors))
  nMods =length(colorLevels)
  nGenes = length(colors);

  gold = substring(opt$MEgold, 3);
  grey = substring(opt$MEgrey, 3);

  # Flag modules whose size is 1
  modSizes = table(colors)
  act = (modSizes>1)

  svds=list()
  kIM = list();
  #printFlush(".coreCalcForAdj:getting svds and kIM");
  if (!opt$densityOnly)
  {            
     svds[[1]] = .getSVDs(datRef, colors);
     kIM[[1]] = .kIM(datRef, colors, calculateAll = opt$calculateCor.kIMall);
  }
   
  if (opt$calculatePermutation | opt$densityOnly)
  {
    svds[[2]] = .getSVDs(datRefP, colors);
    kIM[[2]] = .kIM(datRefP, colors, calculateAll = opt$calculateCor.kIMall);
  } else {
    svds[[2]] = svds[[1]];
    kIM[[2]] = kIM[[1]];
  }

  svds[[3]] = .getSVDs(datTest, colors);
  kIM[[3]] = .kIM(datTest, colors, calculateAll = opt$calculateCor.kIMall);

  proVar=matrix(NA, nMods ,2)

  modGenes = list();
  for (m in 1:nMods)
  {
    modGenes[[m]] = c(1:nGenes)[colors==colorLevels[m]];
    proVar[m, 1] = svds[[2]][[m]]$d[1]/sum(svds[[2]][[m]]$d); 
    proVar[m, 2] = svds[[3]][[m]]$d[1]/sum(svds[[3]][[m]]$d); 
  }

  #printFlush(".coreCalcForAdj:getting corkME and ICOR");
  corkME = rep(NA, nMods);
  corkMEall = rep(NA, nMods);
  corkIM = rep(NA, nMods);
  ICORdat = rep(NA,nMods)
  if (!opt$densityOnly)
  {
     for(m in 1:nMods ) if(act[m])
     {
        nModGenes=modSizes[m];
        corExpr = parse(text=paste(opt$corFnc, "(svds[[1]][[m]]$u,svds[[3]][[m]]$u",
                                                prepComma(opt$corOptions), ")"));
        corkME[m] = abs(eval(corExpr));

        if (opt$calculateCor.kIMall)
        {
          corExpr = parse(text=paste(opt$corFnc, "(kIM[[1]][,m],kIM[[3]][,m] ", 
                                                prepComma(opt$corOptions), ")"));
          corkMEall[m] = eval(corExpr);
        }

        corExpr = parse(text=paste(opt$corFnc, "(kIM[[1]][modGenes[[m]],m],kIM[[3]][modGenes[[m]],m] ", 
                                   prepComma(opt$corOptions), ")"));
        corkIM[m] = eval(corExpr);

        adj1 = datRef[modGenes[[m]], modGenes[[m]]];
        adj2 = datTest[modGenes[[m]], modGenes[[m]]];
        corExpr = parse(text=paste(opt$corFnc, "(c(as.dist(adj1)), c(as.dist(adj2))",
                                   prepComma(opt$corOptions), ")"));
        ICORdat[m] = eval(corExpr);
     }
  }
      
  meanSignAwareKME=matrix(NA, nMods ,2)
  meankIM=matrix(NA, nMods ,2)
  for(m in 1:nMods ) if(act[m])
  {       
     meankIM[m, 1] = mean(kIM[[2]][modGenes[[m]], m], na.rm = TRUE)
     meankIM[m, 2] = mean(kIM[[3]][modGenes[[m]], m], na.rm = TRUE)
     if (opt$densityOnly)
     {
        if (opt$nType==1)
        {
          meanSignAwareKME[m,1]=mean(abs(svds[[2]][[m]]$u),na.rm = TRUE)
          meanSignAwareKME[m,2]=mean(abs(svds[[3]][[m]]$u),na.rm = TRUE)
        } else {
          meanSignAwareKME[m,1]=abs(mean(svds[[2]][[m]]$u,na.rm = TRUE))
          meanSignAwareKME[m,2]=abs(mean(svds[[3]][[m]]$u,na.rm = TRUE))
        }
     } else {
        meanSignAwareKME[m,1]=mean(abs(svds[[2]][[m]]$u),na.rm = TRUE)
        meanSignAwareKME[m,2]=abs(mean(sign(svds[[1]][[m]]$u) * svds[[3]][[m]]$u,na.rm = TRUE))
     }
  }
  
  MeanAdj = matrix(NA, nMods, 2);
  sepMat = array(NA, dim = c(nMods, nMods, 2));
  corCC = rep(NA, nMods);
  corMAR = rep(NA, nMods);
  meanCC = matrix(NA,nMods ,2)
  meanMAR = matrix(NA,nMods ,2)
  for (m in 1:nMods) if (act[m])
  {
    modAdj = datRefP[modGenes[[m]], modGenes[[m]]];
    if (opt$calculateClusterCoeff) ccRefP = .clusterCoeff(modAdj);
    marRefP = .MAR(modAdj);
    meanMAR[m, 1] = mean(marRefP);
    MeanAdj[m,1]=mean(as.dist(modAdj), na.rm = TRUE);

    modAdj = datRef[modGenes[[m]], modGenes[[m]]];
    if (opt$calculateClusterCoeff) ccRef = .clusterCoeff(modAdj);
    marRef = .MAR(modAdj);
  
    modAdj = datTest[modGenes[[m]], modGenes[[m]]];
    if (opt$calculateClusterCoeff) ccTest = .clusterCoeff(modAdj);
    marTest = .MAR(modAdj);
    MeanAdj[m,2] = mean(as.dist(modAdj), na.rm = TRUE);

    if (opt$calculateClusterCoeff)
    {
      meanCC[m, 1] = mean(ccRefP);
      meanCC[m, 2] = mean(ccTest);
      corExpr = parse(text=paste(opt$corFnc, "(ccRef, ccTest ", prepComma(opt$corOptions), ")"));
      corCC[m] = eval(corExpr);
    }
    meanMAR[m, 2] = mean(marTest);
    corExpr = parse(text=paste(opt$corFnc, "(marRef, marTest ", prepComma(opt$corOptions), ")"));
    corMAR[m] = eval(corExpr);

    if ((m > 1) && (colorLevels[m]!=gold))
    {
       for (m2 in 1:(m-1)) if (colorLevels[m2]!=gold)
       {
          interAdj = datRefP[modGenes[[m]], modGenes[[m2]]];
          tmp = mean(interAdj, na.rm = TRUE);
          if (tmp!=0) {
            sepMat[m, m2, 1] = mean(interAdj, na.rm = TRUE)/sqrt(MeanAdj[m, 1] * MeanAdj[m2, 1]);
          } else 
            sepMat[m, m2, 1] = 0;
          sepMat[m2, m, 1] = sepMat[m, m2, 1];
          interAdj = datTest[modGenes[[m]], modGenes[[m2]]];
          tmp = mean(interAdj, na.rm = TRUE);
          if (tmp!=0) {
            sepMat[m, m2, 2] = mean(interAdj, na.rm = TRUE)/sqrt(MeanAdj[m, 2] * MeanAdj[m2, 2]);
          } else
            sepMat[m, m2, 2] = 0;
          sepMat[m2, m, 2] = sepMat[m, m2, 2];
       }
    }
  }
  Separability=matrix(NA, nMods ,2)
  notGold = colorLevels!=gold;
  for(k in 1:2)
     Separability[notGold, k]=1-apply(sepMat[notGold, notGold, k, drop = FALSE], 1, max, na.rm = TRUE)                        
  MeanSignAwareCorDat=matrix(NA,nMods ,2)
  list(modSizes = modSizes, 
       corkIM = corkIM, corkME = corkME, corkMEall = corkMEall, 
       proVar = proVar, 
       meanSignAwareKME = meanSignAwareKME,
       meankIM = meankIM,
       Separability = Separability, MeanSignAwareCorDat = MeanSignAwareCorDat, ICORdat = ICORdat,
       MeanAdj = MeanAdj, meanClusterCoeff = meanCC, meanMAR = meanMAR, corCC = corCC, corMAR = corMAR)
}
