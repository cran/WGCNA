#==========================================================================================================
#
# Utility functions for handling possibly disk-backed blockwise data. 
#
#==========================================================================================================

.getAttributesOrEmptyList = function(object)
{
  att = attributes(object);
  if (is.null(att)) list() else att;
}

newBlockwiseData = function(data, external = FALSE, fileNames = NULL, 
                             doSave = external,
                             recordAttributes = TRUE,
                             metaData = list())
{
  if (length(external)==0) 
     stop("'external' must be logical of length 1.");

  if (!is.null(dim(data)) || !is.list(data))
      stop("'data' must be a list without dimensions.");

  if (recordAttributes)
  {
    attributes = lapply(data, .getAttributesOrEmptyList);
  } else 
    attributes = NULL;

  nBlocks = length(data);

  if (length(metaData) > 0)
  {
     if (length(metaData)!=nBlocks) 
         stop("If 'metaData' are given, it must be a list with one component per component of 'data'.");
  } else {
    metaData = .listRep(list(), nBlocks)
  }

  lengths = sapply(data, length);

  if (doSave && !external)
    warning("newBlockwiseData: Cannot save when 'external' is not TRUE. Data will not be written to disk.")

  if (external)
  {
    if (is.null(fileNames)) stop("When 'external' is TRUE, 'fileNames' must be given.");
  } else
    fileNames = NULL;

  out = list(external = external, data = if (external) list() else data, fileNames = fileNames, 
               lengths = lengths, attributes = attributes, metaData = metaData);
  if (doSave && external)
  {
    if (nBlocks!=length(fileNames)) stop("Length of 'data' and 'fileNames' must be the same.");
    mapply(function(object, f) save(object, file = f), data, fileNames);
  }

  class(out) = "BlockwiseData";
  out;
}


mergeBlockwiseData = function(...)
{
  args = list(...);
  args = args[ sapply(args, length) > 0];
  if (!all(sapply(args, inherits, "BlockwiseData"))) 
    stop("All arguments must be of class 'BlockwiseData'.");
  
  external1 = .checkLogicalConsistency(args, "external");
  .checkListNamesConsistency(lapply(args, getElement, "attributes"), "attributes");
  .checkListNamesConsistency(lapply(args, getElement, "metaData"), "metaData");

  out = list(external = external1, data = do.call(c, lapply(args, .getElement, "data")),
             fileNames = do.call(c, lapply(args, .getElement, "fileNames")),
             lengths = do.call(c, lapply(args, .getElement, "lengths")),
             attributes = do.call(c, lapply(args, .getElement, "attributes")),
             metaData = do.call(c, lapply(args, .getElement, "metaData")));
  class(out) = "BlockwiseData";
  out;
}


# Under normal circumstance arguments external, dist and diag should not be set by the calling fnc, but this
# function can also be used to start a new instance of blockwise data.

addBlockToBlockwiseData = function(bwData, 
               blockData, 
               external = bwData$external,
               blockFile = NULL, 
               doSave = external,
               recordAttributes = !is.null(bwData$attributes),
               metaData = NULL)
{
  badj1 = newBlockwiseData(external = external,
                           data = if (is.null(blockData)) NULL else list(blockData), 
                           fileNames = blockFile, 
                           recordAttributes = recordAttributes,
                           metaData = list(metaData),
                           doSave = doSave)
  mergeBlockwiseData(bwData, badj1);
}

BD.actualFileNames = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure.");
  if (bwData$external) bwData$fileNames else character(0);
}

BD.nBlocks = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure.");
  length(bwData$lengths);
}

  
BD.blockLengths = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise structure.");
  bwData$lengths;
}

BD.getMetaData = function(bwData, blocks = NULL, simplify = TRUE)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise structure.");
  if (is.null(blocks)) blocks = 1:BD.nBlocks(bwData);
  if ( (length(blocks)==0) | any(!is.finite(blocks)))
    stop("'block' must be present and finite.");

  if (any(blocks<1) | (blocks > BD.nBlocks(bwData)))
    stop("All entries in 'block' must be between 1 and ", BD.nBlocks(bwData))

  out = bwData$metaData[blocks];
  if (length(blocks)==1 && simplify)
    out= out[[1]];
  out;
}

.getBDorPlainData = function(data, blocks = NULL, simplify = TRUE)
{
  if (inherits(data, "BlockwiseData")) BD.getData(data, blocks, simplify) else data;
}

BD.getData = function(bwData, blocks = NULL, simplify = TRUE)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure.");

  if (is.null(blocks)) blocks = 1:BD.nBlocks(bwData);
  if ( (length(blocks)==0) | any(!is.finite(blocks)))
    stop("'block' must be present and finite.");

  if (any(blocks<1) | (blocks > BD.nBlocks(bwData))) 
    stop("All entries in 'block' must be between 1 and ", BD.nBlocks(bwData))

  if (bwData$external)
  {
    lengths = BD.blockLengths(bwData);
    out = mapply(.loadObject, bwData$fileNames[blocks], name = 'object', size = lengths[blocks],
                 SIMPLIFY = FALSE);
  } else
    out = bwData$data[blocks];
  if (length(blocks)==1 && simplify)
    out= out[[1]];
  out;
}

BD.checkAndDeleteFiles = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure");
  if (bwData$external)
    .checkAndDelete(bwData$fileNames)
}

.getData = function(x, ...)
{
  if (inherits(x, "BlockwiseData")) return(BD.getData(x, ...));
  x;
}

.setAttr = function(object, name, value)
{
  attr(object, name) = value;
  object;
}

.setAttrFromList = function(object, valueList)
{
  if (length(valueList) > 0) for (i in 1:length(valueList))
      attr(object, names(valueList)[i]) = valueList[[i]];
  object;
}

# A version of getElement that returns NULL if name does not name a valid object
.getElement = function(lst, name)
{
  if (name %in% names(lst)) lst[[name, exact = TRUE]] else NULL
}

.checkLogicalConsistency = function(objects, name)
{
  vals = sapply(objects, getElement, name)
  if (!all(vals) && !all(!vals))
    stop("All arguments must have the same value of '", name, "'.");
  vals[1];
}

.checkListNamesConsistency = function(lst, name)
{
  names = lapply(lst, names);
  if (!all(sapply(names, function(x) isTRUE(all.equal(x, names[[1]])))))
    stop("Not all names agree in ", name);
}

