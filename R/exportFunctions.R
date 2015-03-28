# Functions for exporting networks to various network visualization software

exportNetworkToVisANT = function(
  adjMat,
  file = NULL,
  weighted = TRUE,
  threshold = 0.5,
  maxNConnections = NULL,
  probeToGene = NULL)
{
  adjMat = as.matrix(adjMat)
  adjMat[is.na(adjMat)] = 0;
  nRow = nrow(adjMat);
  checkAdjMat(adjMat, min = -1, max = 1);
  probes = dimnames(adjMat)[[1]]
  if (!is.null(probeToGene))
  {
    probes2genes = match(probes, probeToGene[,1]);
    if (sum(is.na(probes2genes)) > 0) 
      stop("Error translating probe names to gene names: some probe names could not be translated.");
    probes = probeToGene[probes2genes, 2]
  }

  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE);
  colMat = matrix(c(1:nRow), nRow, nRow);

  adjDst = as.dist(adjMat);
  dstRows = as.dist(rowMat);
  dstCols = as.dist(colMat);

  if (is.null(maxNConnections)) maxNConnections = length(adjDst);

  ranks = rank(-abs(adjDst), na.last = TRUE, ties.method = "first")
  edges = abs(adjDst) > threshold & ranks <= maxNConnections
  nEdges = sum(edges)

  visAntData = data.frame (
     from = probes[dstRows[edges]],
     to = probes[dstCols[edges]],
     direction = rep(0, nEdges),
     method = rep("M0039", nEdges),
     weight = if (weighted) adjDst[edges] else rep(1, nEdges)
     );

  if (!is.null(file))
    write.table(visAntData, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  invisible(visAntData);
}

exportNetworkToCytoscape = function(
  adjMat,
  edgeFile = NULL,
  nodeFile = NULL,
  weighted = TRUE,
  threshold = 0.5,
  nodeNames = NULL,
  altNodeNames = NULL,
  nodeAttr = NULL,
  includeColNames = TRUE)
{
  adjMat = as.matrix(adjMat)
  adjMat[is.na(adjMat)] = 0;
  nRow = nrow(adjMat);
  checkAdjMat(adjMat, min = -1, max = 1);
  if (is.null(nodeNames)) nodeNames = dimnames(adjMat)[[1]]
  if (is.null(nodeNames)) 
    stop("Cannot determine node names: nodeNames is NULL and adjMat has no dimnames.")
  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE);
  colMat = matrix(c(1:nRow), nRow, nRow);

  if (!is.null(nodeAttr))
  {
    if (is.null(dim(nodeAttr))) nodeAttr = data.frame(nodeAttribute = nodeAttr);
    nodeAttr = as.data.frame(nodeAttr);
  } else nodeAttr = data.frame(nodeAttribute = rep(NA, ncol(adjMat)));


  adjDst = as.dist(adjMat);
  dstRows = as.dist(rowMat);
  dstCols = as.dist(colMat);

  edges = abs(adjDst) > threshold
  nEdges = sum(edges)
  edgeData = data.frame (
     fromNode = nodeNames[dstRows[edges]],
     toNode = nodeNames[dstCols[edges]],
     weight = if (weighted) adjDst[edges] else rep(1, nEdges),
     direction = rep("undirected", nEdges),
     fromAltName = if (is.null(altNodeNames)) rep("NA", nEdges) else altNodeNames[dstRows[edges]],
     toAltName = if (is.null(altNodeNames)) rep("NA", nEdges) else altNodeNames[dstCols[edges]]
     );

  nodesPresent = rep(FALSE, ncol(adjMat));
  nodesPresent[dstRows[edges]] = TRUE;
  nodesPresent[dstCols[edges]] = TRUE;
  nNodes = sum(nodesPresent);
  nodeData = data.frame (
     nodeName = nodeNames[nodesPresent],
     altName = if (is.null(altNodeNames)) rep("NA", nNodes) else altNodeNames[nodesPresent],
     nodeAttr[nodesPresent, ]
     );

  if (!is.null(edgeFile))
    write.table(edgeData, file = edgeFile, quote = FALSE, row.names = FALSE, col.names = includeColNames,
                sep = "\t");
  
  if (!is.null(nodeFile))
    write.table(nodeData, file = nodeFile, quote = FALSE, row.names = FALSE, col.names = includeColNames,
                sep = "\t");

  list(edgeData = edgeData, nodeData = nodeData);
}
  



