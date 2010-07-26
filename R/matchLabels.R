# Relabel the labels in source such that modules with high overlap with those in reference will have the
# same labels


overlapTable = function(labels1, labels2)
{
  labels1 = as.vector(labels1);
  labels2 = as.vector(labels2);
  levels1 = sort(unique(labels1));
  levels2 = sort(unique(labels2));
  n1 = length(levels1);
  n2 = length(levels2);
  countMat = matrix(0, n1, n2);
  pMat = matrix(0, n1, n2);

  for (m1 in 1:n1)
    for (m2 in 1:n2)
    {
      m1Members = (labels1 == levels1[m1]);
      m2Members = (labels2 == levels2[m2]);
      #print(paste("table for levels", levels1[m1], levels2[m2]));
      #print(table(m1Members, m2Members));
      pMat[m1, m2] = fisher.test(m1Members, m2Members, alternative = "greater")$p.value;
      countMat[m1, m2] = sum(labels1 == levels1[m1] & labels2 == levels2[m2])
    }

 dimnames(pMat) = list(levels1, levels2);
 dimnames(countMat) = list(levels1, levels2);

 pMat[is.na(pMat)] = 1;

 list(countTable = countMat, pTable = pMat)
}
 

matchLabels = function(source, reference, pThreshold = 5e-2 )
{
  if (suppressWarnings(is.na(as.numeric(levels(factor(reference))[1]))))
  {
    origSrc = source;
    origRef = reference;
    stdColorsX = c("grey", standardColors());
    dimSource = dim(source)
    source = match(source, stdColorsX) - 1;
    dim(source) = dimSource;
    reference = match(reference, stdColorsX) -1;
    if (sum(is.na(source)) > 0)
      stop("'source' contains non-numeric labels that do not match standard colors."); 
    if (sum(is.na(reference)) > 0)
      stop("'reference' contains non-numeric labels that do not match standard colors."); 
    colorLabels = TRUE;
  } else
    colorLabels = FALSE;

  source = as.matrix(source);
  if (nrow(source)!=length(reference))
    stop("Number of rows of 'source' must equal the length of 'reference'.");

  result = source;
  refMods = as.numeric(sort(unique(reference)));
  for (col in 1:ncol(source))
  {
    src = source[, col]
    nsource = as.numeric(as.factor(src));
    sourceMods = as.numeric(levels(as.factor(src)));
    newLabels = rep(NA, length(sourceMods));
    tab = overlapTable(src, reference);
    pTab = tab$pTable;
    pOrder = apply(pTab, 2, order);
    bestOrder = order(apply(pTab, 2, min));
    for (rm in 1:length(bestOrder)) if (refMods[bestOrder[rm]]!=0)
    {
      bestInd = 1;
      done = FALSE;
      #printFlush(paste("Looking for best match for reference module ", refMods[bestOrder[rm]]));
      while (!done && bestInd < length(sourceMods))
      {
        bm = pOrder[bestInd, bestOrder[rm]];
        if (sourceMods[bm]!=0)
        {
          bp = pTab[bm, bestOrder[rm]];
          if (bp > pThreshold)
          {
            done = TRUE;
          } else if (is.na(newLabels[bm]))
          {
            #newLabels[bm] = as.numeric(refMods[bestOrder[rm]]);
            #printFlush(paste("Labeling old module ", sourceMods[bm], "as new module",
            #                 refMods[bestOrder[rm]], "with p=", bp));
            newLabels[bm] = refMods[bestOrder[rm]];
            done = TRUE;
          }
        }
        bestInd = bestInd + 1;
      }
    }
    if (is.finite(match(0, sourceMods))) newLabels[match(0, sourceMods)] = 0;
    maxAssd = max(newLabels, refMods, na.rm = TRUE)
    unassdSrcTab = table(source[is.na(newLabels[nsource]), col]);
    unassdRank = rank(-unassdSrcTab);

    newLabels[is.na(newLabels)] = unassdRank + maxAssd;

    result[, col] = newLabels[nsource];
  }

  if (colorLabels) result = labels2colors(result);

  result;
}
    
          
          
        
      


