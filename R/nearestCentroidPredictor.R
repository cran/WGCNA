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
    x, y, 
    xtest = NULL,

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
  doTest = !is.null(xtest)

  if (doTest) 
  {
    xtest = as.matrix(xtest);
    nTestSamples = nrow(xtest);
    if (ncol(x)!=ncol(xtest))
      stop("Number of learning and testing predictors (columns of x, xtest) must equal.");
  } else {
    if (weighSimByPrediction > 0)
       stop("weighting similarity by prediction is not possible when xtest = NULL.");
    nTestSamples = 0;
  }

  numYLevels = sort(unique(y));

  minY = min(y);
  maxY = max(y);

  nSamples = length(y);
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

  # Feature selection:

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
  xSD = apply(x, 2, sd, na.rm = TRUE);
  keep = is.finite(featureSignificance) & (xSD>0);
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
  selectSignif = featureSignificance[select];

  if (scaleFeatureMean)
  {
    if (scaleFeatureVar)
    {
      xSD = apply(xSel, 2, sd, na.rm = TRUE);
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
  if (doTest)
  {
    xtestSel = xtest[, select];
    xtestSel = (xtestSel - matrix(xMean, nTestSamples, nSelect, byrow = TRUE) ) / 
             matrix(xSD, nTestSamples, nSelect, byrow = TRUE);
  } else
    xtestSel = NULL;

  xWeighted = xSel * sampleWeights;

  if (weighSimByPrediction > 0)
  {
    pr = .quickGeneVotingPredictor.CV(xSel, xtestSel, c(1:nSelect))
    dCV = sqrt(colMeans( (pr$CVpredicted - xSel)^2, na.rm = TRUE));
    dTS = sqrt(colMeans( (pr$predictedTest - xtestSel)^2, na.rm = TRUE));
    dTS[dTS==0] = min(dTS[dTS>0]);
    validationWeight = (dCV/dTS)^weighSimByPrediction;
    validationWeight[validationWeight > 1] = 1;
  } else
    validationWeight = rep(1, nSelect);

  nTestSamples = if (doTest) nrow(xtest) else 0;

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
    if (doTest && sum(is.na(xtestSel))>0)
    {
       xtestImp = t(impute.knn(t(xtestSel), k = min(10, nSelect - 1))$data);
    } else
       xtestImp = xtestSel;
  }

  clusterNumbers = rep(1, nLevels);
  sampleModules = list();
  
  # Trivial cluster labels: clusters equal case classes
  for (l in 1:nLevels)
    clusterLabels[[l]] = rep(l, sum(y==levels[l]))

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
      distExpr = spaste( simFnc, "( wcps, wxSel, ", simOptions, ")");
      sample.centroidSim = eval(parse(text = distExpr));

      # Actual prediction: for each sample, calculate distances to centroid profiles
      if (doTest)
      {
        wxtestSel = t(xtestSel) * featureWeight
        distExpr = spaste( simFnc, "( wcps, wxtestSel, ", simOptions, ")");
        testSample.centroidSim = eval(parse(text = distExpr));
      }
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
      wxSel.keepSamples = t(xSel[keepSamples, ]) * featureWeight;
  
      # Back-substitution prediction: 
      
      distExpr = spaste( simFnc, "( wxSel.keepSamples, wxSel, ", simOptions, ")");
      dst = eval(parse(text = distExpr));
      # Test prediction:
      if (doTest)
      {
        wxtestSel = t(xtestSel) * featureWeight
        distExpr = spaste( simFnc, "( wxSel.keepSamples, wxtestSel, ", simOptions, ")");
        dst.test = eval(parse(text = distExpr));
        sample.centroidSim = matrix(0, nClusters, nSamples);
        testSample.centroidSim = matrix(0, nClusters, nTestSamples);
      }
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
   # Save predictions
   predicted = prediction;
   if (doTest)
   {
     testCentroidSimilarities = t(testSample.centroidSim);
     testprediction = cluster2level[apply(testSample.centroidSim, 2, which.max)];
     predictedTest = testprediction;
   }
  #} else 
  #  stop("Prediction for continouos variables is not implemented yet. Sorry!");

  # Reformat output if factors are to be returned
    
  if (returnFactor)
  {
      predicted.out = factor(originalYLevels[[t]][predicted])
      if (doTest) predictedTest.out = factor(originalYLevels[[t]][predictedTest]);
      if (CVfold > 0)
        CVpredicted.out = factor(originalYLevels[[t]][CVpredicted]);
  } else {
    # Turn ordinal predictions into levels of input traits
    predicted.out = originalYLevels[predicted];
    if (doTest) predictedTest.out = originalYLevels[predictedTest];
    if (CVfold > 0)
      CVpredicted.out = originalYLevels[CVpredicted];
  }
  out = list(predicted = predicted.out,
             predictedTest = if (doTest) predictedTest.out else NULL,
             featureSignificance = featureSignificance, 
             selectedFeatures = selectedFeatures,
             centroidProfiles = if (is.null(useQuantile)) centroidProfiles else NULL,
             testSample2centroidSimilarities = if (doTest) testCentroidSimilarities else NULL,
             featureValidationWeights = validationWeight
             )
  if (CVfold > 0)
    out$CVpredicted = CVpredicted.out;

  out;
}


