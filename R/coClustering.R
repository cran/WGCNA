# coClustering and a permutation test for it

coClustering = function(clusters.ref, clusters.test, tupletSize=2, unassignedLabel=0)
{
  overlap = overlapTable(clusters.test, clusters.ref)
  greyRow = rownames(overlap$countTable)==unassignedLabel;
  greyCol = colnames(overlap$countTable)==unassignedLabel;

  refModSizes = table(clusters.ref);

  ccNumer = apply(overlap$countTable[!greyRow, !greyCol, drop = FALSE], 2, choose, tupletSize);
  ccDenom = choose(refModSizes[!greyCol], tupletSize)

  apply(ccNumer, 2, sum)/ccDenom;
}

coClustering.permutationTest = function(clusters.ref, clusters.test, 
                                        tupletSize=2, nPermutations = 100, unassignedLabel=0, 
                                        randomSeed = 12345, verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent);
  if (!is.null(randomSeed)) set.seed(randomSeed);
  observed = coClustering(clusters.ref, clusters.test, tupletSize, unassignedLabel);

  nModules = length(observed)

  permValues = matrix(NA, nPermutations, nModules);

  if (verbose > 0)
    pind = initProgInd(spaste(spaces, "Running permutations: "), " done");
  for (p in 1:nPermutations)
  {
    ctPerm = sample(clusters.test);
    permValues[p, ] = as.numeric(coClustering(clusters.ref, ctPerm, 
                                        tupletSize, unassignedLabel));
    if (verbose > 0) pind = updateProgInd(p/nPermutations, pind);
  }
  if (verbose > 0) printFlush("");
  means = colMeans(permValues);
  sds = apply(permValues, 2, sd, na.rm = TRUE);
  list(observed = observed, Z = (observed-means)/sds, permuted.mean = means, permuted.sd = sds,
       permuted.cc = permValues);
}




