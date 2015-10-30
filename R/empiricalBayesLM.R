#===================================================================================================
#
# Multi-variate empirical Bayes with proper accounting for sigma
#
#===================================================================================================

.colWeightedMeans.x = function(data, weights, na.rm)
{
  nc = ncol(data);
  means = rep(NA, nc);

  for (c in 1:nc)
    means[c] = weighted.mean(data[, c], weights[, c], na.rm = na.rm);

  names(means) = colnames(data);
  means;
}

.weightedScale = function(data, weights)
{
  weightSums = colSums(weights);
  means = .colWeightedMeans.x(data, weights, na.rm = TRUE);
  V1 = colSums(weights);
  V2 = colSums(weights^2);
  centered = data - matrix(means, nrow(data), ncol(data), byrow = TRUE);
  scale = sqrt( colSums(centered^2* weights, na.rm = TRUE) / ( V1 - V2/V1));

  scaled = centered/ matrix(scale,  nrow(data), ncol(data), byrow = TRUE);

  attr(scaled, "scaled:center") = means;
  attr(scaled, "scaled:scale") = scale;
  scaled;
}

.weightedVar = function(x, weights)
{
  V1 = sum(weights);
  V2 = sum(weights^2);
  mean = sum(x * weights)/V1;
  centered = x-mean;
  sum(centered^2 * weights)/(V1-V2/V1);
}

empiricalBayesLM = function(
  data, 
  removedCovariates,
  retainedCovariates = NULL, 
  weights = NULL,
  weightType = c("apriori", "empirical"),
  stopOnSmallWeights = TRUE,
  tol = 1e-4, maxIterations = 1000,
  scaleMeanToSamples = NULL,
  robustPriors = FALSE,
  automaticWeights = c("none", "bicov"),
  aw.maxPOutliers = 0.1)
{

  nSamples = nrow(data);
  designMat = NULL;
  #mean.x = NULL;
  #scale.x = NULL;

  automaticWeights = match.arg(automaticWeights);
  if (automaticWeights=="bicov")
  {
    weightType = "empirical";
    weights = bicovWeights(data, maxPOutliers = aw.maxPOutliers);
  }

  weightType = match.arg(weightType);
  wtype = match(weightType, c("apriori", "empirical"));

  if (!is.null(retainedCovariates))
  {
     if (is.null(dim(retainedCovariates))) 
        retainedCovariates = data.frame(retainedCovariate = retainedCovariates)
     if (any(is.na(retainedCovariates))) stop("All elements of 'retainedCovariates' must be finite.");
     if (nrow(retainedCovariates)!=nSamples)
       stop("Numbers of rows in 'data' and 'retainedCovariates' differ.");
     retainedCovariates = as.data.frame(retainedCovariates);
     mm = model.matrix(~., data = retainedCovariates)[, -1, drop = FALSE];
     colSDs = colSds(mm);
     if (any(colSDs==0))
       stop("Some columns in 'retainedCovariates' have zero variance.");
     designMat = mm;
  }

  if (is.null(removedCovariates))
    stop("'removedCovariates' must be supplied.");
  if (is.null(dim(removedCovariates))) 
      removedCovariates = data.frame(removedCovariate = removedCovariates)
  if (any(is.na(removedCovariates))) stop("All elements of 'removedCovariates' must be finite.");
  if (nrow(removedCovariates)!=nSamples)
    stop("Numbers of rows in 'data' and 'removedCovariates' differ.");
  removedCovariates = as.data.frame(removedCovariates);
  mm = model.matrix(~., data = removedCovariates)[, -1, drop = FALSE];
  colSDs = colSds(mm);
  if (any(colSDs==0))
    stop("Some columns in 'removedCovariates' have zero variance.");
  designMat = cbind(designMat, mm)
  removedColumns = (ncol(designMat)-ncol(mm) + 1):ncol(designMat);

  y.original = as.matrix(data);
  N.original = ncol(y.original);
  if (any(!is.finite(y.original))) 
  {
    warning(immediate. = TRUE,
            "Found missing and non-finite data. These will be removed.");
  }
    
  if (is.null(weights))
    weights = matrix(1, nSamples, ncol(y.original));

  if (any(!is.finite(weights)))
    stop("Given 'weights' contain some infinite or missing entries. All weights must be present and finite.");

  if (any(weights<0))
    stop("Given 'weights' contain negative entries. All weights must be non-negative.");

  originalWeights = weights;

  dimnamesY = dimnames(y.original);

  varY = colVars(y.original, na.rm = TRUE);
  varYMissing = is.na(varY);
  varYZero = varY==0;
  varYZero[is.na(varYZero)] = FALSE;
  keepY = !(varYZero | varYMissing);

  y = y.original[, keepY];
  weights = weights[, keepY];
  yFinite = is.finite(y);
  weights[!yFinite] = 0;
  nSamples.y = colSums(yFinite);

  if (weightType == "apriori")
  { 
    # Check weights
    if (any(weights[yFinite]<1))
    {
      if (stopOnSmallWeights)
      {
         stop("When weights are determined 'apriori', small weights are not allowed.\n",
              "Weights must be at least 1. Use weightType='empirical' if weights were determined from data.")
      } else
         warning(immediate. = TRUE,
                 "Small weights found. This can lead to unreliable fit with weight type 'apriori'.\n",
                 "Proceed with caution.");
    }
  }

  N = ncol(y);
  nc = ncol(designMat);

  # Scale y to mean zero and variance 1. This is needed for the prior to make sense.

  if (is.null(scaleMeanToSamples)) scaleMeanToSamples = c(1:nSamples);

  mean.y.target = .colWeightedMeans.x(y[scaleMeanToSamples, ], weights[scaleMeanToSamples, ], na.rm = TRUE);
  y = .weightedScale(y, weights);
  mean.y = attr(y, "scaled:center");
  mean.y.mat = matrix(mean.y, nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
  scale.y = attr(y, "scaled:scale");
  scale.y.mat = matrix(scale.y, nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
  y[!yFinite] = 0;

  # Get the means of the design matrix with respect to all weight vectors.
  V1 = colSums(weights);
  V2 = colSums(weights^2);
  means.dm = t(designMat) %*% weights / matrix(V1, nrow = nc, ncol = N, byrow = TRUE);

  # Ordinary regression to get starting point for beta

  beta.OLS = matrix(NA, nc, N);
  betaValid = matrix(TRUE, nc, N);
  sigma.OLS = rep(NA, N);
  regressionValid = rep(TRUE, N);
  oldWeights = rep(-1, nSamples);
  for (i in 1:N)
  {
    w1 = weights[, i];
    y1 = y[, i];
    if (any(w1!=oldWeights))
    {
      centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE);
      #dmVar = colSds(centeredDM);
      dmVar.w = apply(centeredDM, 2, .weightedVar, weights = w1);
      #keepDM = dmVar > 0 & dmVar.w > 0;
      keepDM = dmVar.w > 0;
      centeredDM.keep = centeredDM[, keepDM];
      xtx = t(centeredDM.keep) %*% (centeredDM.keep * w1);
      xtxInv = try(solve(xtx), silent = TRUE);
      oldWeights = w1;
    }

    if (!inherits(xtx, "try-error"))
    {
      # The following is really (xtxInv %*% xy1), where xy1 = t(centeredDM) %*% (y[,i]*weights[,i])
      beta.OLS[keepDM, i] = colSums(xtxInv * colSums(centeredDM.keep * y1 * w1));
      betaValid[!keepDM, i] = FALSE;
      y.pred = centeredDM.keep %*% beta.OLS[keepDM, i];
      if (weightType=="apriori")
      {
        # Standard calculation of sigma^2 in weighted refression
        sigma.OLS[i] = sum(w1 * (y1 - y.pred)^2)/(sum(yFinite[, i])-nc-1);
      } else {
        xtxw2 =  t(centeredDM.keep) %*% (centeredDM.keep *w1^2);
        sigma.OLS[i] = sum( w1* (y1-y.pred)^2) / (V1[i] - V2[i]/V1[i] - sum(xtxw2 * xtxInv));
      }
    } else {
      regressionValid[i] = FALSE;
      betaValid[, i] = FALSE;
    }
  }
  if (any(!regressionValid))
     warning(immediate. = TRUE,
             "empiricalBayesLM: OLS regression failed in ", sum(!regressionValid), " variables.");

  # beta.OLS has columns corresponding to variables in data, and rows corresponding to columns in x.

  # Debugging...
  if (FALSE)
  {
    fit = lm(y~., data = data.frame(designMat))
    fit2 = lm(y~., data = as.data.frame(centeredDM));

    max(abs(fit2$coefficients[1,]))

    all.equal(c(fit$coefficients[-1, ]), c(fit2$coefficients[-1, ]))
    all.equal(c(fit$coefficients[-1, ]), c(beta.OLS))
    sigma.fit = apply(fit$residuals, 2, var)*(nSamples-1)/(nSamples-1-nc)
    all.equal(as.numeric(sigma.fit), as.numeric(sigma.OLS))

  }

  # Priors on beta : mean and variance
  if (robustPriors)
  {
    if (is.na(beta.OLS[, regressionValid]))
      stop("Some of OLS coefficients are missing. Please use non-robust priors.");
    prior.means = rowMedians(beta.OLS[, regressionValid], na.rm = TRUE);
    prior.covar = .bicov(t(beta.OLS[, regressionValid]));
  } else {
    prior.means = rowMeans(beta.OLS[, regressionValid], na.rm = TRUE);
    prior.covar = cov(t(beta.OLS[, regressionValid]), use = "complete.obs");
  }
  prior.inverse = solve(prior.covar);

  # Prior on sigma: mean and variance (median and MAD are bad estimators since the distribution is skewed)
  sigma.m = mean(sigma.OLS[regressionValid], na.rm = TRUE);
  sigma.v = var(sigma.OLS[regressionValid], na.rm = TRUE);

  # Turn the sigma mean and variance into the parameters of the inverse gamma distribution
  prior.a = sigma.m^2/sigma.v + 2;
  prior.b = sigma.m * (prior.a-1);

  # Calculate the EB estimates
  beta.EB = beta.OLS;
  sigma.EB = sigma.OLS;
  
  for (i in which(regressionValid))
  {
    # Iterate to solve for EB regression coefficients (betas) and the residual variances (sigma)
    # It appears that this has to be done individually for each variable.

    difference = 1;
    iteration = 1;

    keepDM = betaValid[, i];
    beta.old = as.matrix(beta.OLS[keepDM, i]);
    sigma.old = sigma.OLS[i];
    y1 = y[, i];
    w1 = weights[, i];

    centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE);
    centeredDM.keep = centeredDM[, keepDM];
    xtx = t(centeredDM.keep) %*% (centeredDM.keep * w1);
    xtxInv = solve(xtx);

    if (all(keepDM))
    {
      prior.inverse.keep = prior.inverse;
    } else
      prior.inverse.keep = solve(prior.covar[keepDM, keepDM]);

    while (difference > tol && iteration <= maxIterations)
    {
      y.pred = centeredDM.keep %*% beta.old;
      if (wtype==1)
      {
        # Apriori weights.
        fin1 = yFinite[, i];
        nSamples1 = sum(fin1);
        sigma.new = (sum(w1*(y1-y.pred)^2) + 2*prior.b)/ (nSamples1-nc + 2 * prior.a + 1);
      } else {
        # Empirical weights
        V1 = sum(w1);
        V2 = sum(w1^2);
        xtxw2 =  t(centeredDM.keep) %*% (centeredDM.keep *w1^2);
        sigma.new = (sum( w1* (y1-y.pred)^2) + 2*prior.b) / (V1 - V2/V1 - sum(xtxw2 * xtxInv) + 2*prior.a + 2);
      }

      A = (prior.inverse.keep + xtx/sigma.new)/2;
      A.inv = solve(A);
      #B = as.numeric(t(y1*w1) %*% centeredDM/sigma.new) + as.numeric(prior.inverse %*% prior.means);
      B = colSums(centeredDM.keep * y1 * w1/sigma.new) + colSums(prior.inverse.keep * prior.means[keepDM]);

      #beta.new = A.inv %*% as.matrix(B)/2
      beta.new = colSums(A.inv * B)/2  # ...a different and hopefully faster way of writing the above

      difference = max( abs(sigma.new-sigma.old)/(sigma.new + sigma.old),
                        abs(beta.new-beta.old)/(beta.new + beta.old));

      beta.old = beta.new;
      sigma.old = sigma.new;
      iteration = iteration + 1;
    }
    if (iteration > maxIterations) warning(immediate. = TRUE, 
                                           "Exceeded maximum number of iterations for variable ", i, ".");
    beta.EB[keepDM, i] = beta.old;
    sigma.EB[i] = sigma.old;
  }

  # Put output together. Will return the coefficients for lm and EB-lm, and the residuals with added mean.

  fitAndCoeffs = function(beta, sigma)
  {
      #fitted.removed = fitted = matrix(NA, nSamples, N);
      fitted.removed = fitted = y;
      beta.fin = beta;
      beta.fin[is.na(beta)] = 0;
      for (i in which(regressionValid))
      {
         centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE);
         fitted.removed[, i] = centeredDM[, removedColumns, drop = FALSE] %*% 
                                   beta.fin[removedColumns, i, drop = FALSE];
         fitted[, i] = centeredDM %*% beta.fin[, i, drop = FALSE]
      }
      #browser()
      
      residuals = (y - fitted.removed) * scale.y.mat;
      # Residuals now have weighted column means equal zero.

      meanShift =  matrix(mean.y.target, nSamples, N, byrow = TRUE)

      residuals.all = residualsWithMean.all = fitted.all = matrix(NA, nSamples, N.original);

      residuals.all[, keepY] = residuals;
      residuals.all[, varYZero] = 0;
      residuals.all[is.na(y.original)] = NA;

      residualsWithMean.all[, keepY] = residuals.all[, keepY] + meanShift;
      residualsWithMean.all[is.na(y.original)] = NA;

      fitted.all[, keepY] = fitted * scale.y.mat + meanShift;

      beta.all = beta.all.scaled = matrix(NA, nc+1, N.original);
      sigma.all = sigma.all.scaled = rep(NA, N.original);
      sigma.all.scaled[keepY] = sigma;
      sigma.all[keepY] = sigma * scale.y^2;

      beta.all[-1, keepY] = beta.fin * matrix(scale.y, nrow = nc, ncol = N, byrow = TRUE);
      beta.all.scaled[-1, keepY] = beta.fin;
      beta.all[varYZero] = beta.all.scaled[-1, varYZero] = 0;
      #alpha = mean.y - t(as.matrix(mean.x)) %*% beta * scale.y
      alpha = mean.y - colSums(beta * means.dm, na.rm = TRUE) * scale.y
      beta.all[1, keepY] = alpha;
      beta.all.scaled[1, keepY] = 0;

      dimnames(residualsWithMean.all) = dimnamesY;
      colnames(beta.all) = colnames(beta.all.scaled) = colnames(y.original);
      rownames(beta.all) = rownames(beta.all.scaled) = c("(Intercept)", colnames(designMat));

      list(residuals = residuals.all,
           residualsWithMean = residualsWithMean.all,
           beta = beta.all,
           beta.scaled = beta.all.scaled,
           sigmaSq = sigma.all,
           sigmaSq.scaled = sigma.all.scaled,
           fittedValues = fitted.all);
  }

  fc.OLS = fitAndCoeffs(beta.OLS, sigma.OLS);
  fc.EB = fitAndCoeffs(beta.EB, sigma.EB);

  betaValid.all = matrix(FALSE, nc+1, N.original);
  betaValid.all[-1, keepY] = betaValid;
  betaValid.all[1, keepY] = TRUE;
  dimnames(betaValid.all) = dimnames(fc.OLS$beta);

  list( adjustedData = fc.EB$residualsWithMean,
       residuals = fc.EB$residuals,
       coefficients = fc.EB$beta,
       coefficients.scaled = fc.EB$beta.scaled,
       sigmaSq = fc.EB$sigmaSq,
       sigmaSq.scaled = fc.EB$sigmaSq.scaled,
       fittedValues = fc.EB$fittedValues,

       # OLS results
       adjustedData.OLS = fc.OLS$residualsWithMean,
       residuals.OLS = fc.OLS$residuals,
       coefficients.OLS = fc.OLS$beta,
       coefficients.OLS.scaled = fc.OLS$beta.scaled,
       sigmaSq.OLS = fc.OLS$sigmaSq,
       sigmaSq.OLS.scaled = fc.OLS$sigmaSq.scaled,
       fittedValues.OLS = fc.OLS$fittedValues,


       # Weights used in the model
       weights = originalWeights,

       # indices of valid fits
       dataColumnValid = keepY,
       dataColumnWithZeroVariance = varYZero,
       coefficientValid = betaValid.all);
}


#===================================================================================================
#
# Robust functions
#
#===================================================================================================

.bicov = function(x)
{
  x = as.matrix(x);
  if (any(!is.finite(x))) 
    stop("All entries of 'x' must be finite.");
  nc = ncol(x);
  nr = nrow(x);
  mx = colMedians(x);
  mx.mat = matrix(mx, nr, nc, byrow = TRUE);
  mads = colMads(x, constant = 1);
  mad.mat = matrix(mads, nr, nc, byrow = TRUE);
  u = abs((x - mx.mat)/(9 * mad.mat));
  a = matrix(as.numeric(u<1), nr, nc);

  topMat = a * (x-mx.mat) * (1-u^2)^2;
  top = nr * t(topMat) %*% topMat;
  botVec = colSums(a * (1-u^2) * (1-5*u^2));
  bot = botVec %o% botVec;

  out = top/bot;
  dimnames(out) = list(colnames(x), colnames(x));
  out;
}


# Argument refWeight:
# w = (1-u^2)^2
# u^2 = 1-sqrt(w)
# referenceU = sqrt(1-sqrt(referenceW))

bicovWeights = function(x, pearsonFallback = TRUE, maxPOutliers = 1,
                        outlierReferenceWeight = 0.5625,
                        defaultWeight = 0)
{
  referenceU = sqrt(1-sqrt(outlierReferenceWeight^2));
  dimX = dim(x);
  dimnamesX = dimnames(x);
  x = as.matrix(x);
  nc = ncol(x);
  nr = nrow(x);
  mx = colMedians(x, na.rm = TRUE);
  mx.mat = matrix(mx, nr, nc, byrow = TRUE);
  mads = colMads(x, constant = 1, na.rm = TRUE);
  madZero = mads==0;
  if (any(madZero)) 
  {
     warning(immediate. = TRUE,
             "MAD is zero in some columns of 'x'.");
     if (pearsonFallback)
     {
       sds = colSds(x[, madZero, drop = FALSE], na.rm = TRUE)
       mads[madZero] = sds * qnorm(0.75);
     }
  }
  mad.mat = matrix(mads, nr, nc, byrow = TRUE);
  u = (x - mx.mat)/(9 * mad.mat);
  if (maxPOutliers < 0.5)
  {
    lq = colQuantileC(u, p = maxPOutliers);
    uq = colQuantileC(u, p = 1-maxPOutliers);
    lq[is.na(lq)] = 0;
    uq[is.na(uq)] = 0;

    lq[lq>-referenceU] = -referenceU;
    uq[uq < referenceU] = referenceU;
    lq = abs(lq);
    changeNeg = which(lq>referenceU);
    changePos = which(uq > referenceU);

    for (c in changeNeg)
    {
      neg1 = u[, c] < 0;
      neg1[is.na(neg1)] = FALSE;
      u[neg1, c] = u[neg1, c] * referenceU/lq[c];
    }

    for (c in changePos)
    {
      pos1 = u[, c] > 0;
      pos1[is.na(pos1)] = FALSE;
      u[pos1, c] = u[pos1, c] * referenceU/uq[c];
    }
  }

  a = matrix(as.numeric(abs(u)<1), nr, nc);
  weights = a * (1-u^2)^2;
  weights[is.na(x)] = 0;
  weights[!is.finite(weights)] = defaultWeight;
  dim(weights) = dimX;
  if (!is.null(dimX)) dimnames(weights) = dimnamesX;
  weights;
}

#=====================================================================================================
#
# Diagnostic plot for coefficients
#
#=====================================================================================================




