
# These functions are writen in the framework where for several sets the expression data are a vector of
# lists, with each list having a component "data" in which the actual expression data for the set are
# stored.

#-----------------------------------------------------------------------------------------------
#
# Overall options and settings for the package
#
#-----------------------------------------------------------------------------------------------

.moduleColorOptions = list(MEprefix = "ME")

moduleColor.getMEprefix = function()
{
  .moduleColorOptions$MEprefix;
}

# ===================================================
#The function moduleEigengenes finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "colors".
# The theoretical underpinnings are described in Horvath, Dong, Yip (2005)
# http://www.genetics.ucla.edu/labs/horvath/ModuleConformity/
# This requires the R library impute

moduleEigengenes = function(expr, colors, impute = TRUE, nPC = 1, align = "along average",
                            excludeGrey = FALSE, grey = ifelse(is.numeric(colors),  0, "grey"),
                            subHubs = TRUE, trapErrors = FALSE, 
                            returnValidOnly = trapErrors,
                            softPower = 6, scale = TRUE,
                            verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent);

  if (verbose==1) 
     printFlush(paste(spaces, "moduleEigengenes: Calculating", nlevels(as.factor(colors)), 
                              "module eigengenes in given set."));
  if (is.null(expr))
  {  
    stop("moduleEigengenes: Error: expr is NULL. ");
  }
  if (is.null(colors))
  {  
    print("moduleEigengenes: Error: colors is NULL. ");
    stop()
  }

  if (is.null(dim(expr)) || length(dim(expr))!=2)
    stop("moduleEigengenes: Error: expr must be two-dimensional.");

  if (dim(expr)[2]!=length(colors))
    stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).");

  if (is.factor(colors))
  {
    nl = nlevels(colors);
    nlDrop = nlevels(colors[, drop = TRUE]);
    if (nl > nlDrop)
     stop(paste("Argument 'colors' contains unused levels (empty modules). ", 
                "Use colors[, drop=TRUE] to get rid of them."));
  }

  if (softPower < 0) stop("softPower must be non-negative");

  alignRecognizedValues =  c("", "along average");
  if (!is.element(align, alignRecognizedValues)) {
    printFlush(paste("ModulePrincipalComponents: Error:",
                "parameter align has an unrecognised value:", 
                align, "; Recognized values are ", alignRecognizedValues));
    stop()
  }

  maxVarExplained = 10;
  if (nPC>maxVarExplained)
    warning(paste("Given nPC is too large. Will use value", maxVarExplained));

  nVarExplained = min(nPC, maxVarExplained);
  modlevels=levels(factor(colors))
  if (excludeGrey)
    if (sum(as.character(modlevels)!=as.character(grey))>0) {
      modlevels = modlevels[as.character(modlevels)!=as.character(grey)]
    } else {
      stop(paste("Color levels are empty. Possible reason: the only color is grey",
                 "and grey module is excluded from the calculation."));
    }
  PrinComps = data.frame(matrix(NA,nrow=dim(expr)[[1]], ncol= length(modlevels))) 
  averExpr = data.frame(matrix(NA,nrow=dim(expr)[[1]], ncol= length(modlevels))) 
  varExpl= data.frame(matrix(NA, nrow= nVarExplained, ncol= length(modlevels)))
  validMEs = rep(TRUE, length(modlevels));
  validAEs = rep(FALSE, length(modlevels));
  isPC = rep(TRUE, length(modlevels));
  isHub = rep(FALSE, length(modlevels));
  validColors = colors;
  names(PrinComps)=paste(moduleColor.getMEprefix(), modlevels, sep="")
  names(averExpr)=paste("AE",modlevels,sep="")
  for(i in c(1:length(modlevels)) )
  {
    if (verbose>1) 
      printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", modlevels[i]));
    modulename = modlevels[i]
    restrict1 = as.character(colors)== as.character(modulename)
    if (verbose > 2)
       printFlush(paste(spaces, " ...", sum(restrict1), "genes"));
    datModule = as.matrix(t(expr[, restrict1]));
    n = dim(datModule)[1]; p = dim(datModule)[2];
    pc = try( 
      {
        if (nrow(datModule)>1 && impute)
        {
          seedSaved = FALSE;
          if (exists(".Random.seed")) {
             saved.seed = .Random.seed;
             seedSaved = TRUE;
          }
          if (verbose > 5) printFlush(paste(spaces, " ...imputing missing data"));
          datModule = impute.knn(as.matrix(datModule), k = min(10, nrow(datModule)-1))
          # some versions of impute.knn return a list and we need the data component:
          try( { if (!is.null(datModule$data)) datModule = datModule$data; }, silent = TRUE )
          # The <<- in the next line is extremely important. Using = or <- will create a local variable of
          # the name .Random.seed and will leave the important global .Random.seed untouched.
          if (seedSaved) .Random.seed <<- saved.seed;
        }
        if (verbose > 5) printFlush(paste(spaces, " ...scaling"));
        if (scale) datModule=t(scale(t(datModule)));
        if (verbose > 5) printFlush(paste(spaces, " ...calculating SVD"));
        svd1 = svd(datModule, nu = min(n, p, nPC), nv = min(n, p, nPC));
        # varExpl[,i]= (svd1$d[1:min(n,p,nVarExplained)])^2/sum(svd1$d^2)
        if (verbose > 5) printFlush(paste(spaces, " ...calculating PVE"));
        veMat = cor(svd1$v[, c(1:min(n,p,nVarExplained))], t(datModule), use = "p") 
        varExpl[c(1:min(n,p,nVarExplained)),i]= apply(veMat^2, 1, mean, na.rm = TRUE)
        # this is the first principal component
        svd1$v[,1]
      }, silent = TRUE);
    if (class(pc)=='try-error')
    {
      if ( (!subHubs) && (!trapErrors) ) stop(pc);
      if (subHubs)
      {
        if (verbose>0)
        {
          printFlush(paste(spaces, " ..principal component calculation for module", 
                                   modulename, "failed with the following error:"));
          printFlush(paste(spaces, "     ", pc, spaces,
                           " ..hub genes will be used instead of principal components."));
        }
        isPC[i] = FALSE;
        pc = try( 
        {
          scaledExpr = scale(t(datModule));
          covEx = cov(scaledExpr, use = "p");
          modAdj = abs(covEx)^softPower;
          kIM = (apply(modAdj, 1, sum, na.rm = TRUE))^3;
          if (max(kIM, na.rm = TRUE) > 1) kIM = kIM-1;
          kIM[is.na(kIM)] = 0;
          hub = which.max(kIM)
          alignSign = sign(covEx[, hub]);
          alignSign[is.na(alignSign)] = 0;
          isHub[i] = TRUE;
          pcxMat = scaledExpr * 
                matrix(kIM * alignSign, nrow = nrow(scaledExpr), ncol = ncol(scaledExpr), byrow = TRUE) /
                sum(kIM);
          pcx = apply(pcxMat, 1, sum, na.rm = TRUE);
          varExpl[1, i] = mean(cor(pcx, t(datModule), use = "p")^2, na.rm = TRUE)
          pcx
        }, silent = TRUE);
      }
    }
    
    if (class(pc)=='try-error')
    {
      if (!trapErrors) stop(pc);
      if (verbose>0)
      {
        printFlush(paste(spaces, " ..ME calculation of module", modulename, 
                                 "failed with the following error:"));
        printFlush(paste(spaces, "     ", pc, spaces,
                         " ..the offending module has been removed."));
      }
      warning(paste("Eigengene calculation of module", modulename, 
                    "failed with the following error \n     ", 
                    pc, "The offending module has been removed.\n"));
      validMEs[i] = FALSE; isPC[i] = FALSE; isHub[i] = FALSE;
      validColors[restrict1] = grey;
    } else {
      PrinComps[, i] = pc;
      ae = try( 
      {
        if (isPC[i]) scaledExpr = scale(t(datModule));
        averExpr[, i] = apply(scaledExpr, 1, mean, na.rm = TRUE);
        if (align == "along average")
        {
          if (verbose>4) printFlush(paste(spaces,
                          " .. aligning module eigengene with average expression."))
          if (cor(averExpr[,i], PrinComps[,i], use = "p")<0) PrinComps[,i] = -PrinComps[,i]
        }
        0;
      }, silent = TRUE);
      if (class(ae)=='try-error')
      {
        if (!trapErrors) stop(ae);
        if (verbose>0)
        {
          printFlush(paste(spaces, " ..Average expression calculation of module", modulename,
                                   "failed with the following error:"));
          printFlush(paste(spaces, "     ", ae, spaces,
                           " ..the returned average expression vector will be invalid."));
        }
        warning(paste("Average expression calculation of module", modulename,
                      "failed with the following error \n     ",
                      ae, "The returned average expression vector will be invalid.\n"));
      }
      validAEs[i] = !(class(ae)=='try-error');
    }
  } 
  allOK = (sum(!validMEs)==0)
  if (returnValidOnly && sum(!validMEs)>0) 
  {
    PrinComps = PrinComps[, validMEs]
    averExpr = averExpr[, validMEs];
    varExpl = varExpl[, validMEs];
    validMEs = rep(TRUE, times = ncol(PrinComps));
    isPC = isPC[validMEs];
    isHub = isHub[validMEs];
    validAEs = validAEs[validMEs];
  }
  allPC = (sum(!isPC)==0);
  allAEOK = (sum(!validAEs)==0)
  list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varExpl, nPC = nPC, 
       validMEs = validMEs, validColors = validColors, allOK = allOK, allPC = allPC, isPC = isPC,
       isHub = isHub, validAEs = validAEs, allAEOK = allAEOK)
}

#---------------------------------------------------------------------------------------------
#
# removeGrey
#
#---------------------------------------------------------------------------------------------
# This function removes the grey eigengene from supplied module eigengenes.

removeGreyME = function(MEs, greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
{
  newMEs = MEs;
  if (is.vector(MEs) & mode(MEs)=="list")
  {
    warned = 0;
    newMEs = vector(mode = "list", length = length(MEs));
    for (set in 1:length(MEs))
    {
      if (!is.data.frame(MEs[[set]]$data))
        stop("MEs is a vector list but the list structure is missing the correct 'data' component."); 
      newMEs[[set]] = MEs[[set]];
      if (greyMEName %in% names(MEs[[set]]$data))
      {
         newMEs[[set]]$data = MEs[[set]]$data[, names(MEs[[set]]$data)!=greyMEName];
      } else {
         if (warned==0)
         {
           warning("removeGreyME: The given grey ME name was not found among the names of given MEs.");
           warned = 1;
         }
      }
    }
  } else {
    if (length(dim(MEs))!=2) stop("Argument 'MEs' has incorrect dimensions.")
    MEs = as.data.frame(MEs);
    if (greyMEName %in% names(MEs))
    {
       newMEs = MEs[, names(MEs)!=greyMEName];
    } else {
       warning("removeGreyME: The given grey ME name was not found among the names of given MEs.");
    }
  } 
 
  newMEs;
}
#-------------------------------------------------------------------------------------
#
#  ModulePrincipalComponents
#
#-------------------------------------------------------------------------------------
# Has been superseded by moduleEigengenes above.

# ===================================================
# This function collects garbage

collectGarbage=function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}

#--------------------------------------------------------------------------------------
#
# orderMEs
#
#--------------------------------------------------------------------------------------
#
# performs hierarchical clustering on MEs and returns the order suitable for plotting.

orderMEs = function(MEs, greyLast = TRUE, 
                    greyName = paste(moduleColor.getMEprefix(), "grey", sep=""), 
                    orderBy = 1, order = NULL, 
                    useSets = NULL, verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent);

  if ("eigengenes" %in% names(MEs))
  {
     if (is.null(order))
     {
       if (verbose>0) printFlush(paste(spaces, "orderMEs: order not given, calculating using given set", 
                                          orderBy));
       corPC = cor(MEs$eigengenes, use="p")
       disPC = 1-corPC;
       order = .clustOrder(disPC, greyLast = greyLast, greyName = greyName);
     } 
   
     if (length(order)!=dim(MEs$eigengenes)[2])
       stop("orderMEs: given MEs and order have incompatible dimensions.");
    
     orderedMEs = MEs;
     orderedMEs$eigengenes = as.data.frame(MEs$eigengenes[,order]);
     colnames(orderedMEs$eigengenes) = colnames(MEs$eigengenes)[order];
     if (!is.null(MEs$averageExpr))
     {
       orderedMEs$averageExpr = as.data.frame(MEs$averageExpr[, order])
       colnames(orderedMEs$averageExpr) = colnames(MEs$data)[order];
     }
     if (!is.null(MEs$varExplained))
     {
       orderedMEs$varExplained = as.data.frame(MEs$varExplained[, order])
       colnames(orderedMEs$varExplained) = colnames(MEs$data)[order];
     }
     return(orderedMEs);
  } else {
     check = checkSets(MEs, checkStructure = TRUE, useSets = useSets);
     if (check$structureOK)
     {
        multiSet = TRUE;
     } else {
        multiSet = FALSE;
        MEs = fixDataStructure(MEs);
        useSets = NULL; orderBy = 1;
     }
   
     if (!is.null(useSets)) 
       if (is.na(match(orderBy, useSets))) orderBy = useSets[1];
   
     if (is.null(order))
     {
       if (verbose>0) printFlush(paste(spaces, "orderMEs: order not given, calculating using given set", 
                                          orderBy));
       corPC = cor(MEs[[orderBy]]$data, use="p")
       disPC = 1-corPC;
       order = .clustOrder(disPC, greyLast = greyLast, greyName = greyName);
     } 
   
     if (length(order)!=dim(MEs[[orderBy]]$data)[2])
       stop("orderMEs: given MEs and order have incompatible dimensions.");
    
     nSets = length(MEs);
     orderedMEs = MEs;
     if (is.null(useSets)) useSets = c(1:nSets);
     for (set in useSets) 
     {
       orderedMEs[[set]]$data = as.data.frame(MEs[[set]]$data[,order]);
       colnames(orderedMEs[[set]]$data) = colnames(MEs[[set]]$data)[order];
       if (!is.null(MEs[[set]]$averageExpr))
       {
         orderedMEs[[set]]$averageExpr = as.data.frame(MEs[[set]]$averageExpr[, order])
         colnames(orderedMEs[[set]]$averageExpr) = colnames(MEs[[set]]$data)[order];
       }
       if (!is.null(MEs[[set]]$varExplained))
       {
         orderedMEs[[set]]$varExplained = as.data.frame(MEs[[set]]$varExplained[, order])
         colnames(orderedMEs[[set]]$varExplained) = colnames(MEs[[set]]$data)[order];
       }
     }
     if (multiSet) {
       return(orderedMEs);
     } else {
       return(orderedMEs[[1]]$data);
     }
  }
}

#---------------------------------------------------------------------------------------------
#
# .clustOrder
#
#---------------------------------------------------------------------------------------------

.clustOrder = function(distM, greyLast = TRUE, 
                       greyName = paste(moduleColor.getMEprefix(), "grey", sep=""))
{
  distM = as.matrix(distM);
  distNames = dimnames(distM)[[1]];
  greyInd = match(greyName, distNames);
  if (greyLast && !is.na(greyInd)) 
  {
     clusterMEs = (greyName!=distNames);
     if (sum(clusterMEs)>1)
     {
       h = flashClust(as.dist(distM[clusterMEs, clusterMEs]), method = "average");
       order = h$order;
       if (sum(order>=greyInd)>0) order[order>=greyInd] = order[order>=greyInd]+1;
       order = c(order, greyInd);
     } else if (ncol(distM)>1) 
     {
       if (greyInd==1)
       {
         order = c(2, 1)
       } else order = c(1, 2);
     } else order = 1;
  } else {
     if (length(distM)>1)
     {
       h = flashClust(as.dist(distM), method = "average");
       order = h$order;
     } else order = 1;
  }
  order;

 # print(paste("names:", names(distM), collapse = ", "));
 # print(paste("order:", order, collapse=", "))
}
 
#---------------------------------------------------------------------------------------------
#
# consensusOrderMEs
#
#---------------------------------------------------------------------------------------------
# Orders MEs by the dendrogram of their consensus dissimilarity.

consensusOrderMEs = function(MEs, useAbs = FALSE, useSets = NULL, greyLast = TRUE, 
                             greyName = paste(moduleColor.getMEprefix(), "grey", sep=""), 
                             method = "consensus")
{
  # Debugging code:
  #printFlush("consensusOrderMEs:");
  #size = checkSets(MEs);
  #print(size);
  # end debuging code
  Diss = consensusMEDissimilarity(MEs, useAbs = useAbs, useSets = useSets, method = method);
  order = .clustOrder(Diss, greyLast, greyName);
  #print(order)
  orderMEs(MEs, greyLast = greyLast, greyName = greyName, order = order, useSets = useSets);
} 

#---------------------------------------------------------------------------------------------
#
# consensusMEDissimilarity
#
#---------------------------------------------------------------------------------------------
# This function calcualtes a consensus dissimilarity (i.e., correlation) among sets of MEs (more generally,
# any sets of vectors). 
# CAUTION: when not using absolute value, the minimum similarity will favor the large negative values!

consensusMEDissimilarity = function(MEs, useAbs = FALSE, useSets = NULL, method = "consensus")
{
  methods = c("consensus", "majority");
  m = charmatch(method, methods);
  if (is.na(m))
    stop("Unrecognized method given. Recognized values are", paste(methods, collapse =", "));

  nSets = length(MEs);
  MEDiss = vector(mode="list", length = nSets);
  if (is.null(useSets)) useSets = c(1:nSets);
  for (set in useSets)
  {
    if (useAbs)
    {
        diss = 1-abs(cor(MEs[[set]]$data, use="p"));
    } else
    {
        diss = 1-cor(MEs[[set]]$data, use="p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }

  for (set in useSets)
    if (set==useSets[1])
    {
      ConsDiss = MEDiss[[set]]$Diss;
    } else {
      if (m==1) {
         ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
      } else {
         ConsDiss = ConsDiss + MEDiss[[set]]$Diss;
      }
    }

  if (m==2) ConsDiss = ConsDiss/nSets;

  ConsDiss = as.data.frame(ConsDiss);
  names(ConsDiss) = names(MEs[[useSets[1]]]$data);
  rownames(ConsDiss) = names(MEs[[useSets[1]]]$data);

  ConsDiss;
}

#======================================================================================================
# ColorHandler.R
#======================================================================================================

# A set of global variables and functions that should help handling color names for some 400+ modules.
# A vector called .GlobalStandardColors is defined that holds color names with first few entries 
# being the well-known and -loved colors. The rest is randomly chosen from the color names of R,
# excluding grey colors.

#---------------------------------------------------------------------------------------------------------
#
# .GlobalStandardColors 
#
#---------------------------------------------------------------------------------------------------------
# This code forms a vector of color names in which the first entries are given by BaseColors and the rest
# is "randomly" chosen from the rest of R color names that do not contain "grey" nor "gray".

BaseColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan",
                "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen",
                "darkturquoise", "darkgrey",
                "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", 
                "paleturquoise", "violet", "darkolivegreen", "darkmagenta" );

RColors = colors()[-grep("grey", colors())];
RColors = RColors[-grep("gray", RColors)];
InBase = match(BaseColors, RColors);
ExtraColors = RColors[-c(InBase[!is.na(InBase)])];
nExtras = length(ExtraColors);

# Here is the vector of colors that should be used by all functions:

.GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13*c(1:nExtras) +sin(13*c(1:nExtras))) )] );

standardColors = function(n = NULL)
{
  if (is.null(n)) return(.GlobalStandardColors);
  if ((n>0) && (n<=length(.GlobalStandardColors))) 
  {
    return(.GlobalStandardColors[c(1:n)]);
  } else {
    stop("Invalid number of standard colors requested.");
  }
}

rm(BaseColors, RColors, ExtraColors, nExtras, InBase);

#---------------------------------------------------------------------------------------------------------
#
# normalizeLabels
#
#---------------------------------------------------------------------------------------------------------
# "Normalizes" numerical labels such that the largest group is labeled 1, the next largest 2 etc.
# If KeepZero == TRUE, label zero is preserved.

normalizeLabels = function(labels, keepZero = TRUE)
{
  if (keepZero)
  {
    NonZero = (labels!=0);
  }
  else
  {
    NonZero = rep(TRUE, length(labels));
  }
  f = as.numeric(factor(labels[NonZero]));
  t = table(labels[NonZero]);
  # print(t)
  r = rank(-as.vector(t), ties.method = "first");
  norm_labs = rep(0, times = length(labels));
  norm_labs[NonZero] = r[f];
  norm_labs;
}
  
#---------------------------------------------------------------------------------------------------------
#
# labels2colors
#
#---------------------------------------------------------------------------------------------------------
# This function converts integer numerical labels labels into color names in the order either given by
# colorSeq,
# or (if colorSeq==NULL) by standardColors(). If GreyIsZero == TRUE, labels 0 will be assigned
# the color grey; otherwise presence of labels below 1 will trigger an error.
# dimensions of labels (if present) are preserved.

labels2colors = function(labels, zeroIsGrey = TRUE, colorSeq = NULL)
{
  if (is.null(colorSeq)) colorSeq = standardColors();
  if (zeroIsGrey) minLabel = 0 else minLabel = 1

  if (sum(is.na(labels)) > 0)
    stop("'labels' must not contain NA's");
  if (sum( labels>=minLabel )!= length(labels))
    stop(paste("Input error: something's wrong with labels. Either they are not a numeric vector,",
               "or some values are below", minLabel));

  if (max(labels) > length(colorSeq))
  {
     nRepeats = as.integer((max(labels)-1)/length(colorSeq)) + 1;
     warning(paste("labels2colors: Number of labels exceeds number of avilable colors.", 
                   "Some colors will be repeated", nRepeats, "times."))
     extColorSeq = colorSeq;
     for (rep in 1:nRepeats) 
       extColorSeq = c(extColorSeq, paste(colorSeq, ".", rep, sep=""));
  } else {
     nRepeats = 1;
     extColorSeq = colorSeq;
  }
  colors = rep("grey", length(labels));
  colors[labels!=0] = extColorSeq[labels[labels!=0]];
  if (!is.null(dim(labels)))
  {
    dim(colors) = dim(labels);
  }
  colors;
}

#========================================================================================
#
# MergeCloseModules
#
#========================================================================================

#---------------------------------------------------------------------------------
#
# moduleNumber
#
#---------------------------------------------------------------------------------
# Similar to modulecolor2 above, but returns numbers instead of colors, which is oftentimes more useful.
# 0 means unassigned.
# Return value is a simple vector, not a factor.
# Caution: the module numbers are neither sorted nor sequential, the only guarranteed fact is that grey
# probes are labeled by 0 and all probes belonging to the same module have the same number.

moduleNumber = function(dendro, cutHeight = 0.9, minSize = 50)
{
  Branches = cutree(dendro, h = cutHeight);
  NOnBranches = table(Branches);
  TrueBranch = NOnBranches >= minSize;
  Branches[!TrueBranch[Branches]] = 0;
  
  Branches;
}

#--------------------------------------------------------------------------------------
#
# fixDataStructure
#
#--------------------------------------------------------------------------------------
# Check input data: if they are not a vector of lists, put them into the form of a vector of lists.

fixDataStructure = function(data, verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent);
  if ((class(data)!="list") || (class(data[[1]])!="list"))
  {
    if (verbose>0)
      printFlush(paste(spaces, 
                       "fixDataStructure: data is not a vector of lists: converting it into one."));
    x = data;
    data = vector(mode = "list", length = 1);
    data[[1]] = list(data = x);
    rm(x);
  }
  data;
}

#-------------------------------------------------------------------------------------------
#
# checkSets
#
#-------------------------------------------------------------------------------------------
# Checks sets for consistency and returns some diagnostics.

checkSets = function(data, checkStructure = FALSE, useSets = NULL)
{
  nSets = length(data);
  if (is.null(useSets)) useSets = c(1:nSets);
  if (nSets<=0) stop("No data given.");
  structureOK = TRUE;
  if ((class(data)!="list") || (class(data[[useSets[1]]])!="list"))
  {
    if (checkStructure)
    {
      structureOK = FALSE;
      nGenes = 0; nSamples = 0;
    } else {
      stop("data does not appear to have the correct format. Consider using fixDataStructure",
           "or setting checkStructure = TRUE when calling this function.");
    }
  } else {
    nSamples = vector(length = nSets);
    nGenes = dim(data[[useSets[1]]]$data)[2];
    for (set in useSets) 
    {
      if (nGenes!=dim(data[[set]]$data)[2])
      {
        if (checkStructure) 
        {
           structureOK = FALSE;
        } else {
           stop(paste("Incompatible number of genes in set 1 and", set));
        }
      }
      nSamples[set] = dim(data[[set]]$data)[1];
    }
  }

  list(nSets = nSets, nGenes = nGenes, nSamples = nSamples, structureOK = structureOK);
}


#--------------------------------------------------------------------------------------
#
# multiSetMEs
#
#--------------------------------------------------------------------------------------

multiSetMEs = function(exprData, colors, universalColors = NULL, useSets = NULL, 
                       useGenes = NULL, impute = TRUE, 
                       nPC = 1, align = "along average", excludeGrey = FALSE, 
                       grey = ifelse(is.null(universalColors), ifelse(is.numeric(colors), 0, "grey"), 
                                     ifelse(is.numeric(universalColors), 0, "grey")),
                       subHubs = TRUE,
                       trapErrors = FALSE, 
                       returnValidOnly = trapErrors,
                       softPower = 6,
                       verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);

  nSets = length(exprData);
  setsize = checkSets(exprData, useSets = useSets);
  nGenes = setsize$nGenes;
  nSamples = setsize$nSamples;

  if (verbose>0) printFlush(paste(spaces,"multiSetMEs: Calculating module MEs."));
  MEs = vector(mode="list", length=nSets);
  consValidMEs = NULL;
  if (!is.null(universalColors))
    consValidColors = universalColors;
  if (is.null(useSets)) useSets = c(1:nSets);
  if (is.null(useGenes))
  { 
    for (set in useSets) {
      if (verbose>0) printFlush(paste(spaces,"  Working on set", as.character(set), "...")); 
      if (is.null(universalColors)) {
        setColors = colors[,set];
      } else {
        setColors = universalColors; 
      }
      setMEs = moduleEigengenes(expr = exprData[[set]]$data,
                            color = setColors, impute = impute, nPC = nPC, align = align, 
                            excludeGrey = excludeGrey, grey = grey,
                            trapErrors = trapErrors, subHubs = subHubs,
                            returnValidOnly = FALSE, softPower = softPower,
                            verbose = verbose-1, indent = indent+1);
      if (!is.null(universalColors) && (!setMEs$allOK))
      {
        if (is.null(consValidMEs)) {
          consValidMEs = setMEs$validMEs;
        } else {
          consValidMEs = consValidMEs * setMEs$validMEs;
        }
        consValidColors[setMEs$validColors!=universalColors] = 
            setMEs$validColors[setMEs$validColors!=universalColors]
      }
      MEs[[set]] = setMEs;
      names(MEs[[set]])[names(setMEs)=='eigengenes'] = 'data';
      # Here's what moduleEigengenes returns:
      #
      #  list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varExpl, nPC = nPC, 
      #       validMEs = validMEs, validColors = validColors, allOK = allOK, allPC = allPC, isPC = isPC,
      #       isHub = isHub, validAEs = validAEs, allAEOK = allAEOK)
    }
  } else {
    for (set in useSets) {
      if (verbose>0) printFlush(paste(spaces,"  Working on set", as.character(set), "...")); 
      if (is.null(universalColors)) {
        setColors = colors[useGenes ,set];
      } else {
        setColors = universalColors[useGenes]; 
      }
      setMEs = moduleEigengenes(expr = exprData[[set]]$data[, useGenes],
                            color = setColors, impute = impute, nPC = nPC, align = align, 
                            excludeGrey = excludeGrey, grey = grey,
                            trapErrors = trapErrors, subHubs = subHubs,
                            returnValidOnly = FALSE, softPower = softPower,
                            verbose = verbose-1, indent = indent+1);
      if (!is.null(universalColors) && (!setMEs$allOK))
      {
        if (is.null(consValidMEs)) {
          consValidMEs = setMEs$validMEs;
        } else {
          consValidMEs = consValidMEs * setMEs$validMEs;
        }
        consValidColors[setMEs$validColors!=universalColors[useGenes]] = 
            setMEs$validColors[setMEs$validColors!=universalColors[useGenes]]
      }
      MEs[[set]] = setMEs;
      names(MEs[[set]])[names(setMEs)=='eigengenes'] = 'data';
    }
  }
  if (!is.null(universalColors))
  {
    for (set in 1:nSets)
    {
      if (!is.null(consValidMEs)) MEs[[set]]$validMEs = consValidMEs;
      MEs[[set]]$validColors = consValidColors;
    }
  } 
  for (set in 1:nSets)
  {
    MEs[[set]]$allOK = (sum(!MEs[[set]]$validMEs)==0); 
    if (returnValidOnly)
    {
      valid = (MEs[[set]]$validMEs > 0);
      MEs[[set]]$data = MEs[[set]]$data[, valid];
      MEs[[set]]$averageExpr = MEs[[set]]$averageExpr[, valid];
      MEs[[set]]$varExplained = MEs[[set]]$varExplained[, valid];
      MEs[[set]]$isPC =  MEs[[set]]$isPC[valid];
      MEs[[set]]$allPC = (sum(!MEs[[set]]$isPC)==0)
      MEs[[set]]$isHub =  MEs[[set]]$isHub[valid];
      MEs[[set]]$validAEs =  MEs[[set]]$validAEs[valid];
      MEs[[set]]$allAEOK = (sum(!MEs[[set]]$validAEs)==0)
      MEs[[set]]$validMEs = rep(TRUE, times = ncol(MEs[[set]]$data));
    }
  }
  
  MEs;
}

#---------------------------------------------------------------------------------------------
#
# MergeCloseModules
#
#---------------------------------------------------------------------------------------------
# This function merges modules whose MEs fall on one branch of a hierarchical clustering tree

mergeCloseModules = function(exprData, colors, cutHeight = 0.2, MEs = NULL, 
                             impute = TRUE,
                             useAbs = FALSE, 
                             iterate = TRUE,
                             relabel = FALSE, colorSeq = NULL, 
                             getNewMEs = TRUE,
                             getNewUnassdME = TRUE,
                             useSets = NULL, checkDataFormat = TRUE, 
                             unassdColor = ifelse(is.numeric(colors), 0, "grey"), 
                             trapErrors = FALSE,
                             verbose = 1, indent = 0)
{

  MEsInSingleFrame = FALSE;
  spaces = indentSpaces(indent);

  #numCols = is.numeric(colors);
  #facCols = is.factor(colors);
  #charCols = is.character(colors);

  origColors = colors;

  colors = colors[, drop = TRUE];

  greyMEname = paste(moduleColor.getMEprefix(), unassdColor, sep="");

  if (verbose>0) printFlush(paste(spaces, 
            "mergeCloseModules: Merging modules whose distance is less than", cutHeight));

  if (verbose>3) printFlush(paste(spaces, 
            "  .. will look for grey label", greyMEname));

  if (!checkSets(exprData, checkStructure = TRUE, useSets = useSets)$structureOK)
  {
    if (checkDataFormat)
    {
      exprData = fixDataStructure(exprData);
      MEsInSingleFrame = TRUE;
    } else {
      stop("Given exprData appear to be misformatted.");
    }
  }

  setsize = checkSets(exprData, useSets = useSets);
  nSets = setsize$nSets;
  
  if (!is.null(MEs))
  {
    checkMEs = checkSets(MEs, checkStructure = TRUE, useSets = useSets);
    if (checkMEs$structureOK)
    {
      if (nSets!=checkMEs$nSets)
        stop("Input error: numbers of sets in exprData and MEs differ.")
      for (set in 1:nSets)
      {
        if (checkMEs$nSamples[set]!=setsize$nSamples[set])
            stop(paste("Number of samples in MEs is incompatible with subset length for set", set));
      }
    } else {
      if (MEsInSingleFrame)
      {
        MEs = fixDataStructure(MEs);
        checkMEs = checkSets(MEs);
      } else {
        stop("MEs do not have the appropriate structure (same as exprData). ");
      }
    }
  }

  if (setsize$nGenes!=length(colors))
    stop("Number of genes in exprData is different from the length of original colors. They must equal.");

  if ((cutHeight <0) | (cutHeight>(1+as.integer(useAbs)))) 
    stop(paste("Given cutHeight is out of sensible range between 0 and", 1+as.integer(useAbs) ));

  done = FALSE; iteration = 1;

  MergedColors = colors;
  ok = try( 
  {
    while (!done)
    {
      if (is.null(MEs)) 
      {
        MEs = multiSetMEs(exprData, colors = NULL, universalColors = colors,
                        useSets = useSets, impute = impute,
                        subHubs = TRUE, trapErrors = FALSE, excludeGrey = TRUE, 
                        grey = unassdColor,
                        verbose = verbose-1, indent = indent+1);
        MEs = consensusOrderMEs(MEs, useAbs = useAbs, useSets = useSets, greyLast = FALSE);
        collectGarbage();
      } else if (nlevels(as.factor(colors))!=checkMEs$nGenes)
      {
        if ((iteration==1) & (verbose>0)) printFlush(paste(spaces, "  Number of given module colors", 
                  "does not match number of given MEs => recalculating the MEs."))
        MEs = multiSetMEs(exprData, colors = NULL, universalColors = colors,
                        useSets = useSets, impute = impute,
                        subHubs = TRUE, trapErrors = FALSE, excludeGrey = TRUE,
                        grey = unassdColor,
                        verbose = verbose-1, indent = indent+1);
        MEs = consensusOrderMEs(MEs, useAbs = useAbs, useSets = useSets, greyLast = FALSE);
        collectGarbage();
      }
      if (iteration==1) oldMEs = MEs;
  
      # Check colors for number of distinct colors that are not grey
  
      colLevs = as.character(levels(as.factor(colors)));
      if  ( length(colLevs[colLevs!=as.character(unassdColor)])<2 )
      {
        printFlush(paste(spaces, 
           "mergeCloseModules: less than two proper modules."));
        printFlush(paste(spaces, " ..color levels are",
            paste(colLevs, collapse = ", ")));
        printFlush(paste(spaces, " ..there is nothing to merge."));
        MergedNewColors = colors;
        MergedColors = colors;
        nOldMods = 1; nNewMods = 1;
        oldTree = NULL; Tree = NULL;
        break;
      }
  
      # Cluster the found module eigengenes and merge ones that are too close _in all sets_.
  
      MEDiss = vector(mode="list", length = nSets);
      if (is.null(useSets)) useSets = c(1:nSets);
      for (set in useSets)
      {
        useMEs = c(1:dim(MEs[[set]]$data)[2])[names(MEs[[set]]$data)!=greyMEname];
        if (length(useMEs)<2) 
          stop("Something is wrong: there are two or more proper modules, but less than two proper",
               "eigengenes. Please check that the grey color label and module eigengene label", 
               "are correct.");
        if (useAbs)
        {
            diss = 1-abs(cor(MEs[[set]]$data[, useMEs], use = "p"));
        } else {
            diss = 1-cor(MEs[[set]]$data[, useMEs], use = "p");
        }
        MEDiss[[set]] = list(Diss = diss);
      }
     
      for (set in useSets)
        if (set==useSets[1])
        {
          ConsDiss = MEDiss[[set]]$Diss;
        } else {
          ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
        }
     
      nOldMods = nlevels(as.factor(colors));
      Tree = flashClust(as.dist(ConsDiss), method = "average");
      if (iteration==1) oldTree = Tree;
      TreeBranches = as.factor(moduleNumber(dendro = Tree, cutHeight = cutHeight, minSize = 1));
      UniqueBranches = levels(TreeBranches);
      nBranches = nlevels(TreeBranches)
      NumberOnBranch = table(TreeBranches);
      MergedColors = colors;
      
      # Merge modules on the same branch
  
      for (branch in 1:nBranches) if (NumberOnBranch[branch]>1)
      {
        ModulesOnThisBranch = names(TreeBranches)[TreeBranches==UniqueBranches[branch]];
        ColorsOnThisBranch = substring(ModulesOnThisBranch, 3);
        if (is.numeric(origColors)) ColorsOnThisBranch = as.numeric(ColorsOnThisBranch);
        if (verbose>3) 
           printFlush(paste(spaces, "  Merging original colors", 
                            paste(ColorsOnThisBranch, collapse=", ")));
        for (color in 2:length(ColorsOnThisBranch))
          MergedColors[MergedColors==ColorsOnThisBranch[color]] = ColorsOnThisBranch[1];
      }

      MergedColors = MergedColors[, drop = TRUE];
      
      nNewMods = nlevels(as.factor(MergedColors));
  
      if (nNewMods<nOldMods & iterate)
      {
        colors = MergedColors;
        MEs = NULL;
      } else {
        done = TRUE;
      }
      iteration = iteration+1;
    } 
    if (relabel) 
    {
       RawModuleColors = levels(as.factor(MergedColors));
       # relabel the merged colors to the usual order based on the number of genes in each module
       if (is.null(colorSeq)) 
       {
         if (is.numeric(origColors)) {
           colorSeq = c(1:length(table(origColors)));
         } else {
           nNewColors = length(RawModuleColors);
           colorSeq = labels2colors(c(1:nNewColors))
         }
       }
       
      # nGenesInModule = rep(0, nNewMods);
      # for (mod in 1:nNewMods) nGenesInModule[mod] = sum(MergedColors==RawModuleColors[mod]);
       nGenesInModule = table(MergedColors);
     
       SortedRawModuleColors = RawModuleColors[order(-nGenesInModule)]
    
       # Change the color names to the standard sequence, but leave grey grey 
       # (that's why rank in general does not equal color)
       MergedNewColors = MergedColors;
       if (is.factor(MergedNewColors)) MergedNewColors = as.character(MergedNewColors);
       if (verbose>3) printFlush(paste(spaces, "   Changing original colors:"));
       rank = 0;
       for (color in 1:length(SortedRawModuleColors)) if (SortedRawModuleColors[color]!=unassdColor)
       {
         rank = rank + 1;
         if (verbose>3) printFlush(paste(spaces, "      ", SortedRawModuleColors[color], 
                                    "to ", colorSeq[rank]));
         MergedNewColors[MergedColors==SortedRawModuleColors[color]] = colorSeq[rank];
       }
       if (is.factor(MergedColors)) MergedNewColors = as.factor(MergedNewColors);
    } else {
       MergedNewColors = MergedColors; 
    }
    MergedNewColors = MergedNewColors[, drop = TRUE];

    if (getNewMEs)
    {
      if (nNewMods<nOldMods | relabel | getNewUnassdME)
      {
        if (verbose>0) printFlush(paste(spaces, "  Calculating new MEs..."));
        NewMEs = multiSetMEs(exprData, colors = NULL, universalColors = MergedNewColors,
                             useSets = useSets, impute = impute, subHubs = TRUE, trapErrors = FALSE,
                             excludeGrey = !getNewUnassdME, grey = unassdColor,
                             verbose = verbose-1, indent = indent+1);
        newMEs = consensusOrderMEs(NewMEs, useAbs = useAbs, useSets = useSets, greyLast = TRUE,
                                   greyName = greyMEname);
        MEDiss = vector(mode="list", length = nSets);
        useMEs = c(1:dim(newMEs[[1]]$data)[2])[names(newMEs[[1]]$data)!=greyMEname];
        if (length(useMEs)>1) 
        {
          for (set in useSets)
          {
            if (useAbs)
            {
                diss = 1-abs(cor(newMEs[[set]]$data[, useMEs], use = "p"));
            } else {
                diss = 1-cor(newMEs[[set]]$data[, useMEs], use = "p");
            }
            MEDiss[[set]] = list(Diss = diss);
          }
          for (set in useSets)
            if (set==useSets[1])
            {
              ConsDiss = MEDiss[[set]]$Diss;
            } else {
              ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
            }
          Tree = flashClust(as.dist(ConsDiss), method = "average");
        } else Tree = NULL;
      } else {
        newMEs = MEs;
      }
    } else {
       newMEs = NULL;
    }
    if (MEsInSingleFrame) 
    {
      newMEs = newMEs[[1]]$data;
      oldMEs = oldMEs[[1]]$data;
    }
  }, silent = TRUE);

  if (class(ok)=='try-error')
  {
    if (!trapErrors) stop(ok);
    if (verbose>0)
    {
      printFlush(paste(spaces, "Warning: merging of modules failed with the following error:"));
      printFlush(paste('   ', spaces, ok));
      printFlush(paste(spaces, " --> returning unmerged modules and *no* eigengenes."));
    }
    warning(paste("mergeCloseModules: merging of modules failed with the following error:\n",
                  "    ", ok, " --> returning unmerged modules and *no* eigengenes.\n"));
    list(colors = origColors, allOK = FALSE);
  } else {
    list(colors = MergedNewColors, dendro = Tree, oldDendro = oldTree, cutHeight = cutHeight, 
         oldMEs = oldMEs, newMEs = newMEs, allOK = TRUE);
  }
}

  
