# nearest centroid predictor

# 09: Add weighted measures of sample-centroid similarity.
#     Fix CVfold such that leave one out cross-validation works.

# 08: change the way sample clustering is called. Use a do.call.

# 07: modify bagging and add boosting
#     - revert the predictor to state where it can only take a single trait and a single number of features.
#     - make sure the predictor respects nFeatures=0

# -06:
#    - add new arguments assocCut.hi and assocCut.lo
#    - make nNetworkFeatures equal nFeatures
#    = Bagging and boosting is broken.

# 05: add bagging

# 04: add a self-tuning version

# version 03: try to build a sample network predictor. 

# Cluster the samples in each class separately and identify clusters. It's a bit of a question whether we
# can automate the cluster identification completely. Anyway, use the clusters as additional centroids and
# as prediction use the class of each centroid.


# version 02: add the option to use a quantile of class distances instead of distance from centroid




# Inverse distance between colunms of x and y
# assume that y has few columns and compute the distance matrix between columns of x and columns of y.

.euclideanDist.forNCP = function(x, y, use = 'p') 
{ 
  x = as.matrix(x);
  y = as.matrix(y);

  ny = ncol(y)
  diff = matrix(NA, ncol(x), ny);
  
  for (cy in 1:ny)
    diff[, cy] = apply( (x - y[, cy])^2, 2, sum, na.rm = TRUE);

  -diff;
}

.hclust.cutreeDynamic = function(dst, ...)
{
  args = list(...);
  method = "a";
  if (!is.null(args$method)) 
  {
    method = args$method;
    args$method = NULL;
  }
  if (!is.null(args$cutreeMethod))
  {
    args$method = args$cutreeMethod;
    args$cutreeMethod = NULL;
  }
  tree = hclust(dst, method = method);
  args$dendro = tree;
  args$distM = as.matrix(dst);
  do.call(cutreeDynamic, args);
}

.pam.NCP = function(dst, ...)
{
  cluster::pam(dst, ...)$clustering;
}


.sampleClusters = function(
   data, 
   adjFnc = "cor", adjOptions = "use = 'p'",
   networkPower = 2,
   clusteringFnc = ".hclust.cutreeDynamic",
   clusteringOptions = list(minClusterSize = 4, deepSplit = 2, verbose = 0)
  )
{
  x = try( {
    if (nrow(data) < 4) return(rep(1, nrow(data)));
    if (adjFnc == "dist")
    {
      # printFlush("Using euclidean adjacency...");
      dst = as.matrix(dist(data));
    } else {
      corr = eval(parse(text = paste(adjFnc, "(t(data), ", adjOptions, ")")));
      dst =  1 - ((1+corr)/2)^networkPower;
    }
    nNA = colSums(is.na(dst))
    keepSamples = nNA == 0;
    clusteringOptions$dst = as.dist(dst[keepSamples, keepSamples]);
    labels = rep(0, nrow(data));
    labels[keepSamples] = do.call(clusteringFnc, clusteringOptions);
    if (sum(labels!=0)==0) labels[labels==0] = 1;
  } );
  # if (class(x)=="try-error") browser();
  labels;
}


    #blockwiseModules(t(xx), networkType = "signed", corType = "pearson", TOMType = "none", 
    #                    power = networkPower,
    #                    deepSplit = deepSplit, detectCutHeight = NULL, minModuleSize = minClusterSize, 
    #                    numericLabels = TRUE, mergeCutHeight = 1-maxModuleCor, reassignThreshold = 0,
    #                    pamRespectsDendro = FALSE, minKMEtoJoin = 1, minCoreKME = 0, 
    #                    minKMEtoStay = 0)
  

#=====================================================================================================
#
# Main predictor function.
#
#=====================================================================================================
# best suited to prediction of factors.

# Classification: for each level find nFeatures.eachSide that best distinguish the level from all other
# levels. Actually this doesn't make much sense since I will have to put all the distinguishing sets
# together to form the profiles, so features that may have no relationship to a level will be added if
# there are more than two levels. I can fix that using some roundabout analysis but for now forget it. 

# To do: check that the CVfold validation split makes sense, i.e. none of the bins contains all observations
# of any class.

# Could also add robust standardization

# Output a measure of similarity to class centroids
# Same thing for the gene voting predictor

# prediction for heterogenous cases: sample network in training cases get modules and eigensamples and
# similarly in the controls then use all centroids for classification between k nearest neighbor and
# nearest centroid

# CAUTION: the function standardizes each gene (unless standardization is turned off), so the sample
# networks may be different from what would be expected from the supplied data. 

# Work internally with just numeric entries. Corresponding levels of the response are saved and restored at
# the very end. 

nearestCentroidPredictor = function(
    # Input training and test data
    x, y, xtest,

    # Feature weights and selection criteria
    featureSignificance = NULL,
    assocFnc = "cor", assocOptions = "use = 'p'",
    assocCut.hi = NULL, assocCut.lo = NULL,
    nFeatures.hi = 10, nFeatures.lo = 10,
    weighFeaturesByAssociation = 0,
    scaleFeatureMean = TRUE, scaleFeatureVar = TRUE,

    # Predictor options 
    centroidMethod = c("mean", "eigensample"),
    simFnc = "cor", simOptions = "use = 'p'",
    useQuantile = NULL,
    sampleWeights = NULL,
    weighSimByPrediction = 0,

    # Data pre-processing
    nRemovePCs = 0, removePCsAfterScaling = TRUE,

    # Sample network options
    useSampleNetwork = FALSE,
    adjFnc = "cor", adjOptions = "use = 'p'",
    networkPower = 2,
    clusteringFnc = ".hclust.cutreeDynamic",
    clusteringOptions = list(minClusterSize = 4, deepSplit = 2, verbose = 0),

    # What should be returned
    CVfold = 0, returnFactor = FALSE,

    # General options
    randomSeed = 12345,
    verbose = 2, indent = 0)
{

  # For now we do not support regression

  centroidMethod = match.arg(centroidMethod);

  if (simFnc=="dist")
  {
    if (verbose > 0)
      printFlush(paste("NearestCentroidPredictor: 'dist' is changed to a suitable", 
                       "Euclidean distance function.\n",
                       "   Note: simOptions will be disregarded."));
    simFnc = ".euclideanDist.forNCP";
    simOptions = "use = 'p'"
  }

  # Convert factors to numeric variables.
  
  ySaved = y;
  #if (classify) 
  #{ 
    originalYLevels = sort(unique(y)) ;
    y = as.numeric(as.factor(y));
  #}

  x = as.matrix(x);
  xtest = as.matrix(xtest);

  numYLevels = sort(unique(y));

  minY = min(y);
  maxY = max(y);

  nSamples = length(y);
  nTestSamples = nrow(xtest);
  nVars = ncol(x);

  if (!is.null(assocCut.hi))
  {
    if (is.null(assocCut.lo)) assocCut.lo = -assocCut.hi;
  }

  spaces = indentSpaces(indent);

  if (!is.null(useQuantile)) 
  {
     if ( (useQuantile < 0) | (useQuantile > 1) )
      stop("If 'useQuantile' is given, it must be between 0 and 1.");
  }

  if (is.null(sampleWeights)) sampleWeights = rep(1, nSamples);

  # If cross-validation is requested, change the whole flow and use a recursive call.
  if (CVfold > 0)
  {
    if (CVfold > nSamples )
    {
      printFlush("'CVfold' is larger than number of samples. Will perform leave-one-out cross-validation.");
      CVfold = nSamples;
    }
    ratio = nSamples/CVfold;
    if (floor(ratio)!=ratio)
    {
      smaller = floor(ratio);
      nLarger = nSamples - CVfold * smaller
      binSizes = c(rep(smaller, CVfold-nLarger), rep(smaller +1, nLarger));
    } else
      binSizes = rep(ratio, CVfold);
    if (!is.null(randomSeed))
    {
      if (exists(".Random.seed"))
      {
        saved.seed = .Random.seed;
        seedSaved = TRUE;
      } else
        seedSaved = FALSE;
      set.seed(randomSeed);
    }

    sampleOrder = sample(1:nSamples);
      
    CVpredicted = rep(NA, nSamples);
    CVbin = rep(0, nSamples);

    if (verbose > 0) 
    {
      cat(paste(spaces, "Running cross-validation: "));
      if (verbose==1) pind = initProgInd() else printFlush("");
    }

    if (!is.null(featureSignificance)) 
      printFlush(paste("Warning in nearestCentroidPredictor: \n", 
                       "   cross-validation will be biased if featureSignificance was derived", 
                       "from training data."));

    ind = 1;
    for (cv in 1:CVfold)
    {
      if (verbose > 1) printFlush(paste("..cross validation bin", cv, "of", CVfold));
      end = ind + binSizes[cv] - 1;
      samples = sampleOrder[ind:end];
      CVbin[samples] = cv;
      xCVtrain = x[-samples, , drop = FALSE];
      xCVtest = x[samples, , drop = FALSE];
      yCVtrain = y[-samples];
      yCVtest = y[samples];
      CVsampleWeights = sampleWeights[-samples];
      pr = nearestCentroidPredictor(xCVtrain, yCVtrain, xCVtest,
                        #classify = classify,
                        featureSignificance = featureSignificance, 
                        assocCut.hi = assocCut.hi,
                        assocCut.lo = assocCut.lo,
                        nFeatures.hi = nFeatures.hi,
                        nFeatures.lo = nFeatures.lo,
                        useQuantile = useQuantile, 
                        sampleWeights = CVsampleWeights,

                        nRemovePCs = nRemovePCs,
                        removePCsAfterScaling = removePCsAfterScaling,

                        useSampleNetwork = useSampleNetwork,
                        adjFnc = adjFnc, adjOptions = adjOptions,
                        networkPower = networkPower,
                        clusteringFnc = clusteringFnc,
                        clusteringOptions = clusteringOptions,

                        CVfold = 0, returnFactor = FALSE,
                        randomSeed = randomSeed,
                        centroidMethod = centroidMethod,
                        assocFnc = assocFnc, assocOptions = assocOptions,
                        scaleFeatureMean = scaleFeatureMean,
                        scaleFeatureVar = scaleFeatureVar,
                        simFnc = simFnc, simOptions = simOptions,
                        weighFeaturesByAssociation = weighFeaturesByAssociation,
                        weighSimByPrediction = weighSimByPrediction,
                        verbose = verbose - 2, indent = indent + 1)

      CVpredicted[samples] = pr$predictedTest;
      ind = end + 1;
      if (verbose==1) pind = updateProgInd(cv/CVfold, pind);
    }
    if (verbose==1) printFlush("");
  }

  if (nrow(x)!=length(y))
    stop("Number of observations in x and y must equal.");

  if (ncol(x)!=ncol(xtest))
    stop("Number of learning and testing predictors (columns of x, xtest) must equal.");

  # Remove PCs if requested

  if (!removePCsAfterScaling && nRemovePCs > 0)
  {
    x = removePrincipalComponents(x, nRemovePCs);
    xtest = removePrincipalComponents(xtest, nRemovePCs);
  }

  # Feature selection:

  passedScale = FALSE;
  if (removePCsAfterScaling && nRemovePCs > 0)
  {
    passedScale = TRUE;
    if (scaleFeatureMean)
    {
      if (scaleFeatureVar)
      {
        xSD = sd(x, na.rm = TRUE);
      } else 
        xSD = rep(1, nVars);
      xMean = apply(x, 2, mean, na.rm = TRUE);
    } else {
      if (scaleFeatureVar)
      {
        xSD = sqrt(apply(x^2, 2, sum, na.rm = TRUE)) / pmax(apply(!is.na(x), 2, sum) - 1, rep(1, nVars));
      } else
        xSD = rep(1, nVars);
      xMean = rep(0, nVars);
    }
    x = removePrincipalComponents(scale(x, center = xMean, scale = xSD), n = nRemovePCs);
    xtest = removePrincipalComponents(scale(xtest, center = xMean, scale = xSD), nRemovePCs);
  }

    

  xWeighted = x * sampleWeights;
  yWeighted = y * sampleWeights;

  if (is.null(featureSignificance))
  {
    corEval = parse(text = paste(assocFnc, "(xWeighted, yWeighted, ", assocOptions, ")"));
    featureSignificance = as.vector(eval(corEval));
  } else {
    if (length(featureSignificance)!=nVars)
      stop("Given 'featureSignificance' has incorrect length (must be nFeatures).");
  }
  nGood = nVars;
  nNA = sum(is.na(featureSignificance));
  testCentroidSimilarities = list();

  keep = is.finite(featureSignificance);
  nKeep = sum(keep);
  keepInd = c(1:nVars)[keep];
  order = order(featureSignificance[keep]);
  levels = sort(unique(y));
  nLevels = length(levels);
  if (is.null(assocCut.hi))
  {
    nf = c(nFeatures.hi, nFeatures.lo);
    if (nf[2] > 0) ind1 = c(1:nf[2]) else ind1 = c();
    if (nf[1] > 0) ind2 = c((nKeep-nf[1] + 1):nKeep) else ind2 = c();
    indexSelect = unique(c(ind1, ind2));
    if (length(indexSelect) < 1) 
      stop("No features were selected. At least one of 'nFeatures.hi', 'nFeatures.lo' must be nonzero.");
    indexSelect = indexSelect[indexSelect > 0];
    select = keepInd[order[indexSelect]];
  } else {
    indexSelect = (1:nKeep)[ featureSignificance[keep] >= assocCut.hi |
                               featureSignificance[keep] <= assocCut.lo ]
    if (length(indexSelect)<2) 
       stop(paste("'assocCut.hi'", assocCut.hi, "and assocCut.lo", assocCut.lo, 
                  "are too stringent, less than 3 features were selected.\n",
                  "Please relax the cutoffs.")); 
    select = keepInd[indexSelect];
  }
  if ((length(select) < 3) && (simFnc!='dist')) 
  {
       stop(paste("Less than 3 features were selected. Please either relax", 
                  "the selection criteria of use simFnc = 'dist'."));
  }
  selectedFeatures = select;
  nSelect = length(select);

  xSel = x[, select];
  xtestSel = xtest[, select];
  selectSignif = featureSignificance[select];

  if (!passedScale)
  {
    if (scaleFeatureMean)
    {
      if (scaleFeatureVar)
      {
        xSD = sd(xSel, na.rm = TRUE);
      } else 
        xSD = rep(1, nSelect);
      xMean = apply(xSel, 2, mean, na.rm = TRUE);
    } else {
      if (scaleFeatureVar)
      {
        xSD = sqrt(apply(xSel^2, 2, sum, na.rm = TRUE)) / pmax(apply(!is.na(xSel), 2, sum) - 1, rep(1, nSelect));
      } else
        xSD = rep(1, nSelect);
      xMean = rep(0, nSelect);
    }
    xSel = scale(xSel, center = scaleFeatureMean, scale = scaleFeatureVar);
    xtestSel = (xtestSel - matrix(xMean, nTestSamples, nSelect, byrow = TRUE) ) / 
             matrix(xSD, nTestSamples, nSelect, byrow = TRUE);
  } 
  xWeighted = xSel * sampleWeights;

  if (weighSimByPrediction > 0)
  {
    pr = WGCNA:::.quickGeneVotingPredictor.CV(xSel, xtestSel, c(1:nSelect))
    dCV = sqrt(colMeans( (pr$CVpredicted - xSel)^2, na.rm = TRUE));
    dTS = sqrt(colMeans( (pr$predictedTest - xtestSel)^2, na.rm = TRUE));
    dTS[dTS==0] = min(dTS[dTS>0]);
    validationWeight = (dCV/dTS)^weighSimByPrediction;
    validationWeight[validationWeight > 1] = 1;
  } else
    validationWeight = rep(1, nSelect);

  nTestSamples = nrow(xtest);

  predicted = rep(0, nSamples);
  predictedTest = rep(0, nTestSamples);

  clusterLabels = list();
  clusterNumbers = list();

  if ( (centroidMethod=="eigensample") )
  {
    if (sum(is.na(xSel)) > 0)
    { 
       xImp = t(impute.knn(t(xSel), k = min(10, nSelect - 1))$data);
    } else 
       xImp = xSel;
    if (sum(is.na(xtestSel))>0)
    {
       xtestImp = t(impute.knn(t(xtestSel), k = min(10, nSelect - 1))$data);
    } else
       xtestImp = xtestSel;
  }

  clusterNumbers = rep(1, nLevels);
  sampleModules = list();
  if (useSampleNetwork)
  {
    # Get cluster labels for each class of response in each predictor. This could be made much more
    # efficient for multi-response prediction since the adjacencies only need to be calculated once.
    # The following code assumes that there are at least a few samples of each class in each response,
    # otherwise the network code could fail
    ind = 0;
    for (l in 1:nLevels)
    {
      sampleModules [[l]] = .sampleClusters(xWeighted[y==levels[l], , drop = FALSE],
                                     adjFnc = adjFnc, adjOptions = adjOptions,
                                     networkPower = networkPower, 
                                     clusteringFnc = clusteringFnc,
                                     clusteringOptions = clusteringOptions
                                     );
      if (all(sampleModules[[l]]==0)) sampleModules[[l]] = rep(1, sum(y==levels[l]));
      clusterLabels[[l]] = sampleModules[[l]] + ind;
      clusterLabels[[l]] [sampleModules[[l]] == 0] = 0;
      cls = unique(clusterLabels[[l]]);
      clusterNumbers[l] = length(cls[cls!=0]);
      ind = ind + clusterNumbers[l];
    }
  } else { 
    # Trivial cluster labels: clusters equal case classes
    for (l in 1:nLevels)
      clusterLabels[[l]] = rep(l, sum(y==levels[l]))
    clusterNumbers = rep(1, nLevels);
  }

  nClusters = sum(clusterNumbers);
  centroidSimilarities = array(NA, dim = c(nSamples, nClusters));
  testCentroidSimilarities = array(NA, dim = c(nTestSamples, nClusters));
  #if (classify)
  #{
    cluster2level = rep(c(1:nLevels), clusterNumbers);
    featureWeight = validationWeight; 
    if (is.null(useQuantile))
    {
      # Form centroid profiles for each cluster and class
      centroidProfiles = array(0, dim = c(nSelect, nClusters));
      for (cl in 1:nClusters)
      {
        l = cluster2level[cl];
        clusterSamples = c(1:nSamples)[ y==l ] [ clusterLabels[[l]]==cl ];
        if (centroidMethod=="mean")
        {
          centroidProfiles[, cl] = apply(xSel[clusterSamples, , drop = FALSE], 
                                                    2, mean, na.rm = TRUE);
        } else if (centroidMethod=="eigensample")
        {
          cp = svd(xSel[clusterSamples,], nu = 0, nv = 1)$v[, 1];
          cor = cor(t(xSel[clusterSamples,]), cp);
          if (sum(cor, na.rm = TRUE) < 0) cp = -cp;
          centroidProfiles[, cl] = cp;
        }
      }
      if (weighFeaturesByAssociation > 0)
        featureWeight = featureWeight * sqrt(abs(selectSignif)^weighFeaturesByAssociation);
      # Back-substitution prediction: 
      wcps = centroidProfiles * featureWeight;
      wxSel = t(xSel) * featureWeight;
      wxtestSel = t(xtestSel) * featureWeight
      distExpr = spaste( simFnc, "( wcps, wxSel, ", simOptions, ")");
      sample.centroidSim = eval(parse(text = distExpr));

      # Actual prediction: for each sample, calculate distances to centroid profiles
      distExpr = spaste( simFnc, "( wcps, wxtestSel, ", simOptions, ")");
      testSample.centroidSim = eval(parse(text = distExpr));
    } else {
      labelVector = y;
      for (l in 1:nLevels)
        labelVector[y==l] = clusterLabels[[l]];
      keepSamples = labelVector!=0;
      nKeepSamples = sum(keepSamples);
      keepLabels = labelVector[keepSamples];
      if (weighFeaturesByAssociation > 0)
        featureWeight = featureWeight * sqrt(abs(selectSignif)^weighFeaturesByAssociation);
      wxSel = t(xSel) * featureWeight;
      wxtestSel = t(xtestSel) * featureWeight
      wxSel.keepSamples = t(xSel[keepSamples, ]) * featureWeight;
  
      # Back-substitution prediction: 
      
      distExpr = spaste( simFnc, "( wxSel.keepSamples, wxSel, ", simOptions, ")");
      dst = eval(parse(text = distExpr));
      # Test prediction:
      distExpr = spaste( simFnc, "( wxSel.keepSamples, wxtestSel, ", simOptions, ")");
      dst.test = eval(parse(text = distExpr));
      sample.centroidSim = matrix(0, nClusters, nSamples);
      testSample.centroidSim = matrix(0, nClusters, nTestSamples);
      for (l in 1:nClusters)
      {
         #x = try ( {
         lSamples = c(1:nKeepSamples)[keepLabels==l];
         sample.centroidSim[l, ] = colQuantileC(dst[lSamples, ], 1-useQuantile);
         testSample.centroidSim[l, ] = colQuantileC(dst.test[lSamples, ], 1-useQuantile);
         #} )
         #if (class(x) == 'try-error') browser(text = "zastavka.");
      }
    }
    centroidSimilarities = t(sample.centroidSim);
    prediction = cluster2level[apply(sample.centroidSim, 2, which.max)];
    testCentroidSimilarities = t(testSample.centroidSim);
    testprediction = cluster2level[apply(testSample.centroidSim, 2, which.max)];
    # Save predictions
    predicted = prediction;
    predictedTest = testprediction;
  #} else 
  #  stop("Prediction for continouos variables is not implemented yet. Sorry!");

  # Reformat output if factors are to be returned
    
  if (returnFactor)
  {
      predicted.out = factor(originalYLevels[[t]][predicted])
      predictedTest.out = factor(originalYLevels[[t]][predictedTest]);
      if (CVfold > 0)
        CVpredicted.out = factor(originalYLevels[[t]][CVpredicted]);
  } else {
    # Turn ordinal predictions into levels of input traits
    predicted.out = originalYLevels[predicted];
    predictedTest.out = originalYLevels[predictedTest];
    if (CVfold > 0)
      CVpredicted.out = originalYLevels[CVpredicted];
  }
  out = list(predicted = predicted.out,
             predictedTest = predictedTest.out,
             featureSignificance = featureSignificance, 
             selectedFeatures = selectedFeatures,
             centroidProfiles = if (is.null(useQuantile)) centroidProfiles else NULL,
             testSample2centroidSimilarities = testCentroidSimilarities,
             featureValidationWeights = validationWeight
             )
  if (CVfold > 0)
    out$CVpredicted = CVpredicted.out;

  if (useSampleNetwork)
    out$sampleClusterLabels = clusterLabels;
  out;
}


