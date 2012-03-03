#======================================================================================================
#
# multiData.eigengeneSignificance
#
#======================================================================================================
             
multiData.eigengeneSignificance = function(multiData, multiTrait, moduleLabels,
                        multiEigengenes = NULL, 
                        useModules = NULL,
                        corAndPvalueFnc = corAndPvalue, corOptions = list(),
                        corComponent = "cor", getQvalues = FALSE,
                        setNames = NULL, excludeGrey = TRUE,
                        greyLabel = ifelse(is.numeric(moduleLabels), 0, "grey"))
{
  corAndPvalueFnc = match.fun(corAndPvalueFnc);

  size = checkSets(multiData);
  nSets = size$nSets;
  nGenes = size$nGenes;
  nSamples = size$nSamples;

  if (is.null(multiEigengenes)) {
    multiEigengenes = multiSetMEs(multiData, universalColors = moduleLabels, verbose = 0,
                                  excludeGrey = excludeGrey, grey = greyLabel);
  } else {
    eSize = checkSets(multiEigengenes);
    if (!isTRUE(all.equal(eSize$nSamples, nSamples)))
      stop("Numbers of samples in multiData and multiEigengenes must agree.");
  }

  if (!is.null(useModules))
  {
    keep = substring(colnames(multiEigengenes[[1]]$data), 3) %in% useModules;
    if (sum(keep)==0)
      stop("Incorrectly specified 'useModules': no such module(s).");
    if (any( ! (useModules %in% substring(colnames(multiEigengenes[[1]]$data), 3))))
      stop("Some entries in 'useModules' do not exist in the module labels or eigengenes.");
    for (set in 1:nSets)
      multiEigengenes[[set]]$data = multiEigengenes[[set]]$data[, keep, drop = FALSE]; 
  }

  modLevels = substring(colnames(multiEigengenes[[1]]$data), 3);
  nModules = length(modLevels);

  MES = p = Z = nObs = array(NA, dim = c(nModules, nSets));

  haveZs = FALSE;
  for (set in 1:nSets)
  {
    corOptions$x = multiEigengenes[[set]]$data;
    corOptions$y = multiTrait[[set]]$data;
    cp = do.call(corAndPvalueFnc, args = corOptions);
    corComp = grep(corComponent, names(cp));
    pComp = match("p", names(cp));
    if (is.na(pComp)) pComp = match("p.value", names(cp));
    if (is.na(pComp)) stop("Function `corAndPvalueFnc' did not return a p-value.");
    MES[, set] = cp[[corComp]]
    p[, set] = cp[[pComp]];
    if (!is.null(cp$Z)) { Z[, set] = cp$Z; haveZs = TRUE}
    if (!is.null(cp$nObs))
    {
       nObs[, set] = cp$nObs;
    } else
       nObs[, set] = t(is.na(multiEigengenes[[set]]$data)) %*% (!is.na(multiTrait[[set]]$data));
  }

  if (is.null(setNames))
     setNames = names(multiData);

  if (is.null(setNames))
     setNames = spaste("Set_", c(1:nSets));

  colnames(MES) = colnames(p) = colnames(Z) = colnames(nObs) = setNames;
  rownames(MES) = rownames(p) = rownames(Z) = rownames(nObs) = colnames(multiEigengenes[[1]]$data);

  if (getQvalues)
  {
    q = apply(p, 2, qvalue.restricted);
    dim(q) = dim(p);
    dimnames(q) = dimnames(p);
  } else q = NULL;

  if (!haveZs) Z = NULL;

  list(eigengeneSignificance = MES,
       p.value = p,
       q.value = q,
       Z = Z, 
       nObservations = nObs)
}

#==============================================================================================
#
# nSets
#
#==============================================================================================

nSets = function(multiData, ...) 
{
  size = checkSets(multiData, ...);
  size$nSets;
}
  
             

