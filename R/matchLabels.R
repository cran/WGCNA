# Relabel the labels in source such that modules with high overlap with those in reference will have the
# same labels


overlapTable = function(labels1, labels2, na.rm = TRUE, ignore = NULL,
                        levels1 = NULL, levels2 = NULL)
{
  labels1 = as.vector(labels1);
  labels2 = as.vector(labels2);
  if (na.rm)
  {
    keep = !is.na(labels1) & !is.na(labels2);
    labels1 = labels1[keep];
    labels2 = labels2[keep];
  }
  if (is.null(levels1)) 
  {
    levels1 = sort(unique(labels1));
    levels1 = levels1[!levels1 %in%ignore];
  }
  if (is.null(levels2))
  {
    levels2 = sort(unique(labels2));
    levels2 = levels2[!levels2 %in%ignore];
  }
  n1 = length(levels1);
  n2 = length(levels2);
  countMat = matrix(0, n1, n2);
  pMat = matrix(0, n1, n2);

  for (m1 in 1:n1)
    for (m2 in 1:n2)
    {
      m1Members = (labels1 == levels1[m1]);
      m2Members = (labels2 == levels2[m2]);
      tab = .table2.allLevels(m1Members, m2Members, levels.x = c(FALSE, TRUE), levels.y = c(FALSE, TRUE));
      #print(paste("table for levels", levels1[m1], levels2[m2]));
      #print(table(m1Members, m2Members));
      pMat[m1, m2] = fisher.test(tab, alternative = "greater")$p.value;
      countMat[m1, m2] = sum(labels1 == levels1[m1] & labels2 == levels2[m2])
    }

 dimnames(pMat) = list(levels1, levels2);
 dimnames(countMat) = list(levels1, levels2);

 pMat[is.na(pMat)] = 1;

 list(countTable = countMat, pTable = pMat)
}
 

matchLabels = function(source, reference, pThreshold = 5e-2, na.rm = TRUE,
                       ignoreLabels = if (is.numeric(reference)) 0 else "grey", 
                       extraLabels = if (is.numeric(reference)) c(1:1000) else standardColors())
{

  source = as.matrix(source);
  if (nrow(source)!=length(reference))
    stop("Number of rows of 'source' must equal the length of 'reference'.");

  result = array(NA, dim = dim(source));
  #refMods = as.numeric(sort(unique(reference)));
  #refMods = refMods[!refMods %in% ignoreLabels];
  for (col in 1:ncol(source))
  {
    src = source[, col]
    tab = overlapTable(src, reference, na.rm = na.rm, ignore = ignoreLabels);
    pTab = tab$pTable;
    pOrder = apply(pTab, 2, order);
    bestOrder = order(apply(pTab, 2, min));

    refMods = colnames(pTab);
    if (is.numeric(reference)) refMods = as.numeric(refMods);
    sourceMods = rownames(pTab);
    newLabels = rep(NA, length(sourceMods));
    names(newLabels) = sourceMods;
    for (rm in 1:length(bestOrder)) 
    {
      bestInd = 1;
      done = FALSE;
      #printFlush(paste("Looking for best match for reference module ", refMods[bestOrder[rm]]));
      while (!done && bestInd < length(sourceMods))
      {
        bm = pOrder[bestInd, bestOrder[rm]];
        bp = pTab[bm, bestOrder[rm]];
        if (bp > pThreshold)
        {
            done = TRUE;
        } else if (is.na(newLabels[bm])) {
            #newLabels[bm] = as.numeric(refMods[bestOrder[rm]]);
            #printFlush(paste("Labeling old module ", sourceMods[bm], "as new module",
            #                 refMods[bestOrder[rm]], "with p=", bp));
            newLabels[bm] = refMods[bestOrder[rm]];
            done = TRUE;
        }
        bestInd = bestInd + 1;
      }
    }
    if (length(ignoreLabels) > 0)
    {
      newLabels.ignore = ignoreLabels;
      names(newLabels.ignore) = ignoreLabels;
      newLabels = c(newLabels.ignore, newLabels);
    }

    unassigned = src %in% names(newLabels)[is.na(newLabels)];
    if (any(unassigned))
    {
      unassdSrcTab = table(src[!src %in% names(newLabels)]);
      unassdRank = rank(-unassdSrcTab, ties.method = "first");

      nExtra = sum(is.na(newLabels));
      newLabels[is.na(newLabels)] = extraLabels[ !extraLabels %in% 
                                       c(refMods, ignoreLabels, names(newLabels))] [1:nExtra]; 
    }

    result[, col] = newLabels[match(src, names(newLabels))];
  }

  result;
}
    
          
          
        
      


