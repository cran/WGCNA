#========================================================================================================
#
# Convenience functions for manipulating multiData structures
#
#========================================================================================================

# Note: many of these function would be simpler to use if I used some sort of class/method technique to keep
# track of the class of each object internally. For example, I could then write a generic function "subset" that
# would work consistently on lists and multiData objects. Similarly, multiData2list would simply become a
# method of as.list, and as.list would be safe to use both on lists and on multiData objects. 

mtd.subset = function(multiData, rowIndex = NULL, colIndex = NULL, permissive = FALSE, drop = FALSE)
{
  size = checkSets(multiData, checkStructure = permissive);
  if (!size$structureOK && !is.null(colIndex))
    warning(immdiate. = TRUE,
            paste("mtd.subset: applying column selection on data sets that do not have\n",
                  " the same number of columns. This is treacherous territory; proceed with caution."));
  if (is.null(colIndex)) colIndex.1 = c(1:size$nGenes) else colIndex.1 = colIndex;
  if (is.null(rowIndex)) rowIndex = lapply(size$nSamples, function(n) {c(1:n)})
  if (length(rowIndex)!=size$nSets) 
    stop("If given, 'rowIndex' must be a list of the same length as 'multiData'.");
  out = list();
  for (set in 1:size$nSets)
  {
    if (permissive)
      if (is.null(colIndex)) colIndex.1 = c(1:ncol(multiData[[set]]$data)) else colIndex.1 = colIndex;
    out[[set]] = list(data = multiData[[set]]$data[rowIndex[[set]], colIndex.1, drop = drop]);
  }
  names(out) = names(multiData);
  out;
}

multiData2list = function(multiData)
{
  lapply(multiData, `[[`, 'data');
}

list2multiData = function(data)
{
  out = list();
  for (set in 1:length(data))
    out[[set]] = list(data = data[[set]]);
  names(out) = names(data);
  out;
}

mtd.colnames = function(multiData)
{
  colnames(multiData[[1]]$data);
}

mtd.apply = function(multiData, FUN, ..., mdaSimplify = FALSE, mdaCopyNonData = FALSE)
{
  if (!isMultiData(multiData, strict = FALSE))
    stop("Supplied 'multiData' is not a valid multiData structure.");

  if (mdaSimplify && mdaCopyNonData) 
    stop("Non-data copying is not compatible with simplification.");

  nSets = length(multiData);
  if (mdaCopyNonData) out = multiData else out = vector(mode = "list", length = nSets);
  FUN = match.fun(FUN);
  for (set in 1:nSets)
    out[[set]]$data = FUN(multiData[[set]]$data, ...)

  if (mdaSimplify) return (mtd.simplify(out));

  out;
}

mtd.applyToSubset = function(multiData, FUN, ..., 
                             mdaRowIndex = NULL, mdaColIndex = NULL, 
                             mdaSimplify = FALSE, mdaCopyNonData = FALSE)
{
  size = checkSets(multiData);
  if (mdaSimplify && mdaCopyNonData)
    stop("Non-data copying is not compatible with simplification.");

  if (mdaCopyNonData) res = multiData else res = vector(mode = "list", length = size$nSets);

  doSelection = FALSE;
  if (!is.null(mdaColIndex))
  {
    doSelection = TRUE
    if (any(mdaColIndex < 0 | mdaColIndex > size$nGenes)) 
      stop("Some of the indices in 'mdaColIndex' are out of range.");
  } else {
    mdaColIndex = c(1:size$nGenes);
  }

  if (!is.null(mdaRowIndex))
  {
    if (!is.list(mdaRowIndex))
       stop("mdaRowIndex must be a list, with one component per set.");
    if (length(mdaRowIndex)!=size$nSets)
       stop("Number of components in 'mdaRowIndex' must equal number of sets.");
    doSelection = TRUE
  } else {
    mdaRowIndex = lapply(size$nSamples, function(n) { c(1:n) });
  }
    
  fun = match.fun(FUN) 
  for (set in 1:size$nSets)
    res[[set]] = list(data = fun( 
             if (doSelection) multiData[[set]] $ data[mdaRowIndex[[set]], mdaColIndex, drop = FALSE] else
                              multiData[[set]] $ data, ...));
  if (mdaSimplify) return (mtd.simplify(res));

  return(res);
}

mtd.simplify = function(multiData)
{
  len = length(multiData[[1]]$data);
  dim = dim(multiData[[1]]$data);
  simplifiable = TRUE;
  nSets = length(multiData);
  for (set in 1:nSets)
  {
    if (len!=length(multiData[[set]]$data)) simplifiable = FALSE;
    if (!isTRUE(all.equal( dim, dim(multiData[[set]]$data)))) simplifiable = FALSE;
  }
  if (simplifiable)
  {
    if (is.null(dim)) {
       innerDim = len;
       innerNames = names(multiData[[1]]$data);
       if (is.null(innerNames)) innerNames = spaste("X", c(1:len));
    } else {
       innerDim = dim;
       innerNames = dimnames(multiData[[1]]$data);
       if (is.null(innerNames)) 
         innerNames = lapply(innerDim, function(x) {spaste("X", 1:x)})
       nullIN = sapply(innerNames, is.null);
       if (any(nullIN))
         innerNames[nullIN] = lapply(innerDim[nullIN], function(x) {spaste("X", 1:x)})
    }
    setNames = names(multiData);
    if (is.null(setNames)) setNames = spaste("Set_", 1:nSets);
    mtd.s = matrix(NA, prod(innerDim), nSets);
    for (set in 1:nSets)
      mtd.s[, set] = as.vector(multiData[[set]]$data);

    dim(mtd.s) = c(innerDim, nSets);
    if (!is.null(innerNames))
      dimnames(mtd.s) = c (if (is.list(innerNames)) innerNames else list(innerNames), list(setNames));
    return(mtd.s);
  }
  return(multiData);
}

isMultiData = function(x, strict = TRUE)
{
  if (strict) {
     !inherits(try(checkSets(x), silent = TRUE), 'try-error');
  } else {
    hasData = sapply(x, function(l) { "data" %in% names(l) });
    all(hasData)
  }
}

mtd.mapply = function(FUN, ..., MoreArgs = NULL, mdmaSimplify = FALSE, mdma.doCollectGarbage = FALSE,
                            mdma.argIsMultiData = NULL)
{
  dots = list(...);
  if (length(dots)==0) 
    stop("No arguments were specified. Please type ?mtd.mapply to see the help page.");
  dotLengths = sapply(dots, length);
  if (any(dotLengths!=dotLengths[1]))
    stop(spaste("All arguments to vectorize over must have the same length.\n", 
                "Scalar arguments should be put into the 'MoreArgs' argument.\n",
                "Note: lengths of '...' arguments are: ", paste(dotLengths, collapse = ", ")));
  nArgs = length(dots);
  res = list();
  if (is.null(mdma.argIsMultiData)) mdma.argIsMultiData = sapply(dots, isMultiData, strict = FALSE);

  FUN = match.fun(FUN);
  nSets = dotLengths[1];
  for (set in 1:nSets)
  {
    localArgs = list();
    for (arg in 1:nArgs)
      localArgs[[arg]] = if (mdma.argIsMultiData[arg]) dots[[arg]] [[set]] $ data else dots[[arg]] [[set]];
    names(localArgs) = names(dots);
    res[[set]] = list(data = do.call(FUN, c(localArgs, MoreArgs)));
    if (mdma.doCollectGarbage) collectGarbage();
  }

  names(res) = names(dots[[1]]);

  if (mdmaSimplify)
    return(mtd.simplify(res));

  return(res);
}


mtd.rbindSelf = function(multiData)
{
  size = checkSets(multiData);
  out = NULL;
  colnames = mtd.colnames(multiData);
  for (set in 1:nSets)
  {
    if (!is.null(colnames(multiData[[set]]$data)) && 
        !isTRUE(all.equal(colnames, colnames(multiData[[set]]$data))) )
          colnames(multiData[[set]]$data) = colnames;
    out = rbind(out, multiData[[set]]$data);
  }
  out;
}

mtd.setAttr = function(multiData, attribute, valueList)
{
  size = checkSets(multiData);
  ind = 1;
  for (set in 1:size$nSets)
  {
    attr(multiData[[set]]$data, attribute) = valueList[[ind]];
    ind = ind + 1;
    if (ind > length(valueList)) ind = 1;
  }
  multiData
}

mtd.setColnames = function(multiData, colnames)
{
  size = checkSets(multiData);
  for (set in 1:size$nSets)
    colnames(multiData[[set]]$data) = colnames
  multiData
}


