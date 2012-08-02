## changes to 1.0.R
# 1. add an argument, specify those covariates that should always be included in all bags.
# 2. change scaling mode
# 3. change nFeaturesInBag default


.forwardSelection = function( x, y, xtest, nCandidateCovariates,  candidateCorFnc, candidateCorOptions, classify, NmandatoryCovariates) 
{  
  indxMiss = apply(is.na(x), 2, sum)>0
  indxVar = apply(x, 2, var, na.rm=T)==0
  x = x[, (!indxMiss & !indxVar), drop=F]
  xtest = xtest[, (!indxMiss & !indxVar), drop=F]

  if (NmandatoryCovariates>0) 
  {
   NmandatoryCovariates = NmandatoryCovariates - sum(is.element(1:NmandatoryCovariates, c(which(indxMiss), which(indxVar))))
  }
  
  candidateCorOptions$x = x
  candidateCorOptions$y = y
  absGS = abs(do.call(candidateCorFnc, candidateCorOptions))
  ## nCandidateCovariates could be smaller than indicated due to missing data in x.
  nCandidateCovariates = min(nCandidateCovariates, ncol(x))

  if (NmandatoryCovariates>0)
  {
    indx = c(1:NmandatoryCovariates, which(rank(-absGS[-(1:NmandatoryCovariates)], ties.method="f")<=(nCandidateCovariates-NmandatoryCovariates))+NmandatoryCovariates)
  }else{
    indx = which(rank(-absGS, ties.method="f")<=nCandidateCovariates)
  }
  x = x[, indx, drop=F]
  xtest = xtest[, indx, drop=F]
  absGS = absGS[indx]

  datSelectedAsCandidates = colnames(x) 
 
  geneMax = which(absGS==max(absGS, na.rm=T))

  if (classify) {
    eval1 = parse(text=paste("glm(y~", ifelse(NmandatoryCovariates>0, paste(colnames(x)[1:NmandatoryCovariates], collapse="+"), colnames(x)[geneMax]), ", data=as.data.frame(x), family=binomial(link='logit'))"))
    lmInit = eval(eval1)
    lmUpper = glm(y~., data=as.data.frame(x), family=binomial(link="logit"))
  } else {
    eval1 = parse(text= paste("lm(y~", ifelse(NmandatoryCovariates>0, paste(colnames(x)[1:NmandatoryCovariates], collapse="+"), colnames(x)[geneMax]), ", data=as.data.frame(x))"))
    lmInit = eval(eval1)
    lmUpper = lm(y~., data=as.data.frame(x))
  }
  
  model = stepAIC(lmInit, 
		scope = list(upper = lmUpper), 
		direction="forward", 
		trace=F)

  datSelectedByForwardRegression = rownames(summary(model)$coef)[-1]
  datCoefOfForwardRegression = summary(model)$coef[-1,1]

  ## outHat is piHat in binary, and yHat in quantitative
  outHat = predict(model, as.data.frame(xtest), type="response")

  out=list(outHat, datSelectedAsCandidates, datSelectedByForwardRegression, datCoefOfForwardRegression)
  out
}



################################################
## main user level function

randomGLMpredictor = function(
  x, y, xtest = NULL, 
  classify = TRUE,
  nBags = 100,
  replace = TRUE,
  nObsInBag = if (replace) nrow(x) else as.integer(0.632 * nrow(x)),
  nFeaturesInBag = ceiling(ifelse(ncol(x)<=10, ncol(x), 
		ifelse(ncol(x)<=300, (0.68-0.0016*ncol(x))*ncol(x), ncol(x)/5))),
  nCandidateCovariates=50,
  candidateCorFnc= cor,
  candidateCorOptions = list(method = "pearson", use="p"),
  mandatoryCovariates = NULL,
  randomSeed = 12345,
  verbose =1)
{

  ySaved = y;
  xSaved = x;
  
  require("MASS")
  require("WGCNA")

  if (classify)
  {
    originalYLevels = sort(unique(y));
    if (length(originalYLevels)>2) stop("Error: the predictor is not suitable for discrete outcomes with more than 2 levels.")
    y = as.numeric(as.factor(y))-1;
    numYLevels = sort(unique(y));
    minY = min(y, na.rm = TRUE);
    maxY = max(y, na.rm = TRUE);
  }else{
    if (!is.numeric(y)) stop("Error: A quantitative outcome should be numeric.")
  }

  x = as.matrix(x);
  x = matrix(as.numeric(x), nrow(x), ncol(x))

  if (length(y)!=nrow(x)) {stop("x and y don't have the same number of observations.")}

  nSamples = length(y);
  nVars = ncol(x);

  nFeaturesInBag = ifelse(ncol(x)<nFeaturesInBag,{print("Warning: ncol(x)<nFeaturesInBag"); ncol(x)}, nFeaturesInBag)
  nCandidateCovariates = ifelse (nCandidateCovariates>nFeaturesInBag, {print("Warning: nCandidateCovariates>nFeaturesInBag. Use nCandidateCovariates=nFeaturesInBag"); nFeaturesInBag}, nCandidateCovariates)

  mandatory = !is.null(mandatoryCovariates)
  if (mandatory  & nCandidateCovariates<=length(mandatoryCovariates))
  {
    stop("Error: number of mandatoryCovariates >= nCandidateCovariates")
  }

  ## rename gene names
  colnames(x) = paste("gene", 1:nVars, sep="")

  ## scale x
  xSD = apply(x, 2, sd, na.rm=T)
  xMean = apply(x, 2, mean, na.rm=T)
  validFeatures =xSD>0
  x = scale(x)
  x[, !validFeatures] = 0


  doTest = !is.null(xtest);
  if (doTest)
  { 
    xtestSaved = xtest;

    xtest = as.matrix(xtest);
    xtest = matrix(as.numeric(xtest), nrow(xtest), ncol(xtest))

    if (ncol(x)!=ncol(xtest))
      stop("Number of learning and testing predictors (columns of x, xtest) must equal.");

    nTestSamples = nrow(xtest);
    colnames(xtest) = colnames(x)

    ## scale xtest

    xtest = (xtest - matrix(xMean, nTestSamples, nVars, byrow=T))/matrix(xSD, nTestSamples, nVars, byrow=T)
    xtest[, !validFeatures] = 0
    
    predictedTestMat =  matrix(NA, nTestSamples, nBags);
    
  } 

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

  predictedMat = matrix(NA, nSamples, nBags);

  datSelectedByForwardRegression = datSelectedAsCandidates = datCoefOfForwardRegression = matrix(0, nBags, nVars)
  bagObsIndx = matrix(NA, nBags, nObsInBag)

  for (p in 1:nBags)
  {
    bagSamples = sample(nSamples, nObsInBag, replace = replace);
    bagObsIndx[p,] = bagSamples
    oob = c(1:nSamples)[-unique(bagSamples)];
    nOOB = length(oob);

    if (mandatory) 
    {
      features = c(mandatoryCovariates, sample((1:nVars)[-mandatoryCovariates], nFeaturesInBag - length(mandatoryCovariates)))
    }else{
      features = sample(1:nVars, nFeaturesInBag);}

    xBag = x[bagSamples, features];
    yBag = y[bagSamples];

    if (doTest)
    {
      xTestBag = rbind(x[oob, features], xtest[, features]);
    } else {
      xTestBag = x[oob, features];
    }

    pr = .forwardSelection(xBag, yBag, xTestBag, nCandidateCovariates = nCandidateCovariates, candidateCorFnc = candidateCorFnc, candidateCorOptions = candidateCorOptions, classify=classify, NmandatoryCovariates = length(mandatoryCovariates));

    predictedMat[oob, p] = pr[[1]][1:nOOB];
    if (doTest) {
      predictedTestMat[ , p] = pr[[1]][(nOOB+1):(nOOB + nTestSamples)];
    }
    datSelectedAsCandidates[p, as.numeric(substring(pr[[2]], 5))] = 1
    indxModel = as.numeric(substring(pr[[3]], 5))
    datSelectedByForwardRegression[p, indxModel] = 1
    datCoefOfForwardRegression[p, indxModel] = pr[[4]]
  }

  if (verbose==1) printFlush("");

  timesSelectedByForwardRegression = apply(datSelectedByForwardRegression, 2, sum)

  rownames(bagObsIndx) = rownames(datSelectedAsCandidates) = rownames(datSelectedByForwardRegression) =   rownames(datCoefOfForwardRegression) = paste("bag",1:nBags,sep="")

  if (!is.null(colnames(xSaved))) {
    colnames(datSelectedAsCandidates) = 
    colnames(datSelectedByForwardRegression) =   
    colnames(datCoefOfForwardRegression) = 
    names(timesSelectedByForwardRegression)  = colnames(xSaved)
  }

  predictedOOB.cont1 = apply(predictedMat, 1, mean, na.rm = TRUE);
  predictedOOB.cont1[apply(is.na(predictedMat), 1, sum)== nBags] = NA;
  if (!is.null(rownames(xSaved))) {
    names(predictedOOB.cont1) = rownames(xSaved)
  }
  out = list(   predictedOOB.cont = predictedOOB.cont1,
		datSelectedAsCandidates = datSelectedAsCandidates,
		datSelectedByForwardRegression = datSelectedByForwardRegression,
		datCoefOfForwardRegression = datCoefOfForwardRegression,
		bagObsIndx = bagObsIndx,
		timesSelectedByForwardRegression = timesSelectedByForwardRegression)
  
  if (doTest) {
    predictedTest.cont1 = apply(predictedTestMat, 1, mean, na.rm = TRUE);
    predictedTest.cont1[apply(is.na(predictedTestMat), 1, sum)== nBags] = NA;
    if (!is.null(rownames(xtestSaved))) {
       names(predictedTest.cont1) = rownames(xtestSaved)
    }

    out$predictedTest.cont = predictedTest.cont1
  }


  if (classify)
  {
    predictedOOB = round(predictedOOB.cont1)
    predictedOOB = originalYLevels[predictedOOB+1]

    predictedOOB.cont = cbind(1-predictedOOB.cont1, predictedOOB.cont1)
    colnames(predictedOOB.cont) = as.character(originalYLevels)

    out$predictedOOB = predictedOOB
    out$predictedOOB.cont = predictedOOB.cont
	
    if (doTest) {
      predictedTest = round(predictedTest.cont1)
      predictedTest = originalYLevels[predictedTest+1]

      predictedTest.cont = cbind(1-predictedTest.cont1, predictedTest.cont1)
      colnames(predictedTest.cont) = as.character(originalYLevels)

      out$predictedTest.cont = predictedTest.cont;
      out$predictedTest = predictedTest;
    }
  } 


  out
}








