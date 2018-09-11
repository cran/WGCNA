#========================================================================================================
#
# Convenience functions for manipulating multiData structures
#
#========================================================================================================

# Note: many of these function would be simpler to use if I used some sort of class/method technique to keep
# track of the class of each object internally. For example, I could then write a generic function "subset" that
# would work consistently on lists and multiData objects. Similarly, multiData2list would simply become a
# method of as.list, and as.list would be safe to use both on lists and on multiData objects. 

mtd.subset = function(multiData, rowIndex = NULL, colIndex = NULL, invert = FALSE, 
                      permissive = FALSE, drop = FALSE)
                      
{
  if (length(multiData)==0) return(NULL);
  size = checkSets(multiData, checkStructure = permissive);
  if (!size$structureOK && !is.null(colIndex))
    warning(immediate. = TRUE,
            paste("mtd.subset: applying column selection on data sets that do not have\n",
                  " the same number of columns. This is treacherous territory; proceed with caution."));
  if (is.null(colIndex)) { if (!invert) colIndex.1 = c(1:size$nGenes)} else colIndex.1 = colIndex;
  if (is.null(rowIndex)) rowIndex = lapply(size$nSamples, function(n) if (invert) numeric(0) else c(1:n))
  if (length(rowIndex)!=size$nSets) 
    stop("If given, 'rowIndex' must be a list of the same length as 'multiData'.");
  out = list();
  if (invert)
  {
     rowIndexAll = lapply(size$nSamples, function(n) c(1:n))
     if (any(sapply(rowIndex, function(i) any(i<0))))
       stop("Negative 'rowIndex' indices cannot be used with 'invert=TRUE'.");
     rowIndex = mapply(setdiff, rowIndexAll, rowIndex, SIMPLIFY = FALSE);
  }
  for (set in 1:size$nSets)
  {
    if (permissive)
      if (is.null(colIndex) && !invert) colIndex.1 = c(1:ncol(multiData[[set]]$data)) else colIndex.1 = colIndex;
    if (is.character(colIndex.1))
    { 
      colIndex.1 = match(colIndex.1, colnames(multiData[[set]]$data));
      n1 = length(colIndex.1)
      if (any(is.na(colIndex.1)))
        stop("Cannot match the following entries in 'colIndex' to column names in set ", set, ":\n",
             paste( colIndex[is.na(colIndex.1)] [1:min(n1, 5)], collapse = ", "),
             if (n1>5) ", ... [output truncated]" else "");
    }
    if (invert)
    {
      if (any(colIndex.1<0)) 
        stop("Negative indices cannot be used with 'invert=TRUE'.");
      colIndex.2 = setdiff(1:ncol(multiData[[set]]$data), colIndex.1);
    } else 
      colIndex.2 = colIndex.1;

    out[[set]] = list(data = multiData[[set]]$data[rowIndex[[set]], colIndex.2, drop = drop]);
  }
  names(out) = names(multiData);
  out;
}

multiData2list = function(multiData)
{
  lapply(multiData, getElement, 'data');
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
  if (length(multiData)==0) return(NULL);
  colnames(multiData[[1]]$data);
}

.calculateIndicator = function(nSets, mdaExistingResults, mdaUpdateIndex)
{
  if (length(mdaUpdateIndex)==0) mdaUpdateIndex = NULL;
  calculate = rep(TRUE, nSets);
  if (!is.null(mdaExistingResults))
  {
    nSets.existing = length(mdaExistingResults);
    if (nSets.existing>nSets)
      stop("Number of sets in 'mdaExistingResults' is higher than the number of sets in 'multiData'.\n",
           "  Please supply a valid 'mdaExistingResults' or NULL to recalculate all results.");
    if (nSets.existing==0)
      stop("Number of sets in 'mdaExistingResults' is zero.\n",
           "  Please supply a valid 'mdaExistingResults' or NULL to recalculate all results.");
    if (is.null(mdaUpdateIndex))
    {
      calculate[1:length(mdaExistingResults)] = FALSE;
    } else {
      if (any(! mdaUpdateIndex %in% c(1:nSets)))
        stop("All entries in 'mdaUpdateIndex' must be between 1 and the number of sets in 'multiData'.");
      calculateIndex = sort(unique(c(mdaUpdateIndex,
                                      if (nSets.existing<nSets) c((nSets.existing+1):nSets) else NULL)));
      calculate[ c(1:nSets)[-calculateIndex] ] = FALSE;
    }
  }

  calculate;
}

  

mtd.apply = function(
    # What to do
    multiData, FUN, ...,

    # Pre-existing results and update options
    mdaExistingResults = NULL, mdaUpdateIndex = NULL,
    mdaCopyNonData = FALSE,

    # Output formatting options
    mdaSimplify = FALSE,
    returnList = FALSE,

    # Internal behaviour options
    mdaVerbose = 0, mdaIndent = 0
)
{
  if (length(multiData)==0) return(NULL);
  printSpaces = indentSpaces(mdaIndent);

  if (!isMultiData(multiData, strict = FALSE))
    stop("Supplied 'multiData' is not a valid multiData structure.");

  if (mdaSimplify && mdaCopyNonData) 
    stop("Non-data copying is not compatible with simplification.");

  nSets = length(multiData);
  if (mdaCopyNonData) out = multiData else out = vector(mode = "list", length = nSets);

  calculate = .calculateIndicator(nSets, mdaExistingResults, mdaUpdateIndex);

  FUN = match.fun(FUN);
  for (set in 1:nSets) 
  {
    if (calculate[set])
    {
      if (mdaVerbose > 0)
        printFlush(spaste(printSpaces, "mtd.apply: working on set ", set)); 
      out[[set]] = list(data = FUN(multiData[[set]]$data, ...))
    } else
      out[set] = mdaExistingResults[set];
  }

  names(out) = names(multiData);

  if (mdaSimplify) 
  {
    if (mdaVerbose > 0)
      printFlush(spaste(printSpaces, "mtd.apply: attempting to simplify...")); 
    return (mtd.simplify(out));
  } else if (returnList) {
    return (multiData2list(out));
  }

  out;
}

mtd.applyToSubset = function(
    # What to do
    multiData, FUN, ...,

    # Which rows and cols to keep
    mdaRowIndex = NULL, mdaColIndex = NULL,

    # Pre-existing results and update options
    mdaExistingResults = NULL, mdaUpdateIndex = NULL,
    mdaCopyNonData = FALSE,

    # Output formatting options
    mdaSimplify = FALSE,
    returnList = FALSE,

    # Internal behaviour options
    mdaVerbose = 0, mdaIndent = 0
)
{
  if (length(multiData)==0) return(NULL);
  printSpaces = indentSpaces(mdaIndent);

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
    
  calculate = .calculateIndicator(nSets, mdaExistingResults, mdaUpdateIndex);

  fun = match.fun(FUN) 
  for (set in 1:size$nSets)
  {
    if (calculate[set])
    {
       if (mdaVerbose > 0)
         printFlush(spaste(printSpaces, "mtd.applyToSubset: working on set ", set));
       res[[set]] = list(data = fun( 
                if (doSelection) multiData[[set]] $ data[mdaRowIndex[[set]], mdaColIndex, drop = FALSE] else
                                 multiData[[set]] $ data, ...));
    } else
       res[set] = mdaExistingResults[set];
  }

  names(res) = names(multiData);

  if (mdaSimplify) 
  {
    if (mdaVerbose > 0)
      printFlush(spaste(printSpaces, "mtd.applyToSubset: attempting to simplify..."));
    return (mtd.simplify(res));
  } else if (returnList) {
    return (multiData2list(res));
  }

  return(res);
}

mtd.simplify = function(multiData)
{
  if (length(multiData)==0) return(NULL);
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

mtd.mapply = function(
  # What to do
  FUN, ..., MoreArgs = NULL,

  # How to interpret the input
  mdma.argIsMultiData = NULL,

  # Copy previously known results?
  mdmaExistingResults = NULL, mdmaUpdateIndex = NULL,

  # How to format output
  mdmaSimplify = FALSE,
  returnList = FALSE,

  # Options controlling internal behaviour
  mdma.doCollectGarbage = FALSE,
  mdmaVerbose = 0, mdmaIndent = 0)

{
  printSpaces = indentSpaces(mdmaIndent);
  dots = list(...);
  if (length(dots)==0) 
    stop("No arguments were specified. Please type ?mtd.mapply to see the help page.");
  dotLengths = sapply(dots, length);
  if (any(dotLengths!=dotLengths[1]))
  {
    tmp = data.frame(name = names(dots), length = dotLengths);
    rownames(tmp) = NULL;
    on.exit(print(tmp));
    stop(spaste("All arguments to vectorize over must have the same length.\n", 
                "Scalar arguments should be put into the 'MoreArgs' argument.\n",
                "Note: lengths of '...' arguments are: "));
  }
  nArgs = length(dots);
  res = list();
  if (is.null(mdma.argIsMultiData)) mdma.argIsMultiData = sapply(dots, isMultiData, strict = FALSE);

  nSets = dotLengths[1];

  calculate = .calculateIndicator(nSets, mdmaExistingResults, mdmaUpdateIndex);

  FUN = match.fun(FUN);
  for (set in 1:nSets)
  {
    if (calculate[set])
    {
      if (mdmaVerbose > 0)
        printFlush(spaste(printSpaces, "mtd.mapply: working on set ", set));

      localArgs = list();
      for (arg in 1:nArgs)
        localArgs[[arg]] = if (mdma.argIsMultiData[arg]) dots[[arg]] [[set]] $ data else dots[[arg]] [[set]];
      names(localArgs) = names(dots);
      res[[set]] = list(data = do.call(FUN, c(localArgs, MoreArgs)));
      if (mdma.doCollectGarbage) collectGarbage();
    } else
      res[set] = mdmaExistingResults[set];
  }

  names(res) = names(dots[[1]]);

  if (mdmaSimplify)
  {
    if (mdmaVerbose > 0)
      printFlush(spaste(printSpaces, "mtd.mapply: attempting to simplify..."));
    return (mtd.simplify(res));
  } else if (returnList) {
    return (multiData2list(res));
  }

  return(res);
}


mtd.rbindSelf = function(multiData)
{
  if (length(multiData)==0) return(NULL);
  size = checkSets(multiData);
  out = NULL;
  colnames = mtd.colnames(multiData);
  for (set in 1:size$nSets)
  {
    if (!is.null(colnames(multiData[[set]]$data)) && 
        !isTRUE(all.equal(colnames, colnames(multiData[[set]]$data))) )
    {
       warning("mtd.rbindSelf: 'colnames' of the first set and set ", set, " do not agree.");
       colnames(multiData[[set]]$data) = colnames;
    }
  }
  do.call(rbind, multiData2list(multiData));
}

mtd.setAttr = function(multiData, attribute, valueList)
{
  if (length(multiData)==0) return(NULL);
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
  if (length(multiData)==0) return(NULL);
  size = checkSets(multiData);
  for (set in 1:size$nSets)
    colnames(multiData[[set]]$data) = colnames
  multiData
}


multiData = function(...)
{
  list2multiData(list(...));
}  


