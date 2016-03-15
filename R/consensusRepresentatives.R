#============================================================================================================
#
# mdx version of collapseRows
#
#============================================================================================================


.cr.absMaxMean <- function(x, robust)
{
   if (robust)
   {
     colMedians(abs(x), na.rm = TRUE)
   } else
     colMeans(abs(x), na.rm = TRUE)
}

.cr.absMinMean <- function(x, robust)
{
   if (robust)
   {
     -colMedians(abs(x), na.rm = TRUE)
   } else
     -colMeans(abs(x), na.rm = TRUE)
}


.cr.MaxMean <- function(x, robust)
{
   if (robust)
   {
     colMedians(x, na.rm = TRUE)
   } else
     colMeans(x, na.rm = TRUE)
}

.cr.MinMean <- function(x, robust)
{
   if (robust)
   {
     -colMedians(x, na.rm = TRUE)
   } else
     -colMeans(x, na.rm = TRUE)
}

.cr.maxVariance <- function(x, robust)
{
   if (robust) 
   {
     colMads(x, na.rm = TRUE)
   } else 
     colSds(x, na.rm = TRUE)

}

.checkConsistencyOfGroupAndColID = function(mdx, colID, group)
{
  colID  = as.character(colID)
  group = as.character(group)
  if (length(colID)!=length(group))
     stop("'group' and 'colID' must have the same length.")

  if (any(duplicated(colID)))
     stop("'colID' contains duplicate entries.")

  rnDat = mtd.colnames(mdx);

  if ( sum(is.na(colID))>0 )
      warning(spaste("The argument colID contains missing data. It is recommend you choose non-missing,\n",
              "unique values for colID, e.g. character strings."))

   if ( sum(group=="",na.rm=TRUE)>0 ){
      warning(paste("group contains blanks. It is strongly recommended that you remove",
                    "these rows before calling the function.\n",
                    "   But for your convenience, the collapseRow function will remove these rows"));
      group[group==""]=NA
   }

   if ( sum(is.na(group))>0 ){
       warning(paste("The argument group contains missing data. It is strongly recommended\n",
              "   that you remove these rows before calling the function. Or redefine group\n",
           "   so that it has no missing data. But for convenience, we remove these data."))
   }   

  if ((is.null(rnDat))&(checkSets(mdx)$nGenes==length(colID)))
  {
    write("Warning: mdx does not have column names. Using 'colID' as column names.","")
    rnDat = colID;
    mdx = mtd.setColnames(mdx, colID);
  }
  if (is.null(rnDat))
     stop("'mdx' does not have row names and \n",
          "length of 'colID' is not the same as # variables in mdx.");

  keepProbes = rep(TRUE, checkSets(mdx)$nGenes);

  if (sum(is.element(rnDat,colID))!=length(colID)){
      write("Warning: row names of input data and probes not identical...","")
      write("... Attempting to proceed anyway. Check results carefully.","")
      keepProbes = is.element(colID, rnDat);
      colID = colID[keepProbes]
      mdx = mtd.subset(mdx, , colID);
      group = group[colID]
  }

  restCols = (group!="" & !is.na(group))
  if (any(!restCols))
  {
    keepProbes[keepProbes] = restCols;
    mdx = mtd.subset(mdx, , restCols);
    group = group[restCols]
    colID = colID[restCols]
    rnDat = rnDat[restCols]
  }

  list(mdx = mdx, group = group, colID = colID, keepProbes = keepProbes);
}



selectFewestConsensusMissing <- function(mdx, colID, group, 
                                 minProportionPresent = 1,
                                 consensusQuantile = 0,
                                 verbose = 0, ...)

{
## For each gene, select the gene with the fewest missing probes, and return the results.
#   If there is a tie, keep all probes involved in the tie.
#   The main part of this function is run only if omitGroups=TRUE
   
   otherArgs = list(...)
   nVars = checkSets(mdx)$nGenes;
   nSamples = checkSets(mdx)$nSamples;
   nSets = length(mdx);

   if ((!"checkConsistency" %in% names(otherArgs)) || otherArgs$checkConsistency)
   {
      cd = .checkConsistencyOfGroupAndColID(mdx, colID, group);
      mdx = cd$mdx;
      group = cd$group;
      colID = cd$colID;
      keep = cd$keepProbes;
   } else
      keep = rep(TRUE, nVars);

   # First, return datET if there is no missing data, otherwise run the function
   if (sum(mtd.apply(mdx, function(x) sum(is.na(x)), mdaSimplify = TRUE))==0) 
     return(rep(TRUE, nVars));
   
   # Set up the variables.
   names(group)     = colID
   probes              = mtd.colnames(mdx);
   genes               = group[probes]
   keepGenes           = rep(TRUE, nVars);
   tGenes              = table(genes)
   checkGenes          = sort(names(tGenes)[tGenes>1])
   presentData         = as.matrix(mtd.apply(mdx, function(x) colSums(is.finite(x)), mdaSimplify = TRUE));

   presentFrac = presentData/matrix(nSamples, nVars, nSets, byrow = TRUE);
   consensusPresentFrac = .consensusCalculation(setTomMat = presentFrac, useMean = FALSE,
                            setWeightMat = NULL,
                            consensusQuantile = consensusQuantile)$consensus;
   
   # Omit all probes with at least omitFrac genes missing
   #keep = consensusPresentFrac > omitFraction
   minProportionPresent = as.numeric(minProportionPresent);

   # Omit relevant genes and return results
   if (minProportionPresent > 0)
   {
      if (verbose) pind = initProgInd();
      for (gi in 1:length(checkGenes))
      {
         g = checkGenes[gi];
         gn            = which(genes==g)
         keepGenes[gn] = (consensusPresentFrac[gn] >= minProportionPresent * max(consensusPresentFrac[gn]))
         if (verbose) pind = updateProgInd(gi/length(checkGenes), pind);
      }
      if (verbose) printFlush("");
   }

   keep[keep] = keepGenes;
   return (keep);
}

# ----------------- Main Function ------------------- #

consensusRepresentatives = function(mdx, 
                            group, colID, 
                            consensusQuantile = 0,
                            method = "MaxMean", 
                            useGroupHubs = TRUE,
                            calibration = c("none", "full quantile"),
                            selectionStatisticFnc = NULL, 
                            connectivityPower=1, 
                            minProportionPresent=1,
                            getRepresentativeData = TRUE,
                            statisticFncArguments = list(),
                            adjacencyArguments = list(),
                            verbose = 2, indent = 0)

# Change in methodFunction: if the methodFunction picks a single representative, it should return it in
# attribute "selectedRepresentative".
# minProportionPresent now gives the fraction of the maximum of present values that will still be included.
# minProportionPresent=1 corresponds to minProportionPresent=TRUE in original collapseRows.

# In connectivity-based collapsing, use simple connectivity, do not normalize. This way the connectivities
# retain a larger spread which should prevent quantile normalization from making big changes and potentially
# suprious changes. 

{

   if (!is.null(dim(mdx)))
   {
     warning("consensusRepresentatives: wrapping matrix-like input into a mdx structure.");
     mdx = multiData(mdx);
   }
   spaces = indentSpaces(indent);
   nSamples= checkSets(mdx)$nSamples
   nSets = length(mdx);

   colnames.in = mtd.colnames(mdx);

   calibration = match.arg(calibration);
   
    ## Test to make sure the variables are the right length.
    #     if not, fix it if possible, or stop.

   cd = .checkConsistencyOfGroupAndColID(mdx, colID, group);
   colID = cd$colID;
   group = cd$group;
   mdx = cd$mdx;
   keepVars = cd$keepProbes

   rnDat = mtd.colnames(mdx);

      
## For each gene, select the gene with the fewest missing probes (if minProportionPresent==TRUE)
##  Also, remove all probes with more than 90% missing data
   
   if (verbose > 0)
     printFlush(spaste(spaces, "..selecting variables with lowest numbers of missing data.."));
   keep = selectFewestConsensusMissing(mdx, colID, group, minProportionPresent, 
                                 consensusQuantile = consensusQuantile, verbose = verbose -1)
   mdx = mtd.subset(mdx, , keep);
   keepVars[keepVars] = keep;

   group = group[keep];
   colID = colID[keep];

   rnDat = mtd.colnames(mdx);
   
##   If method="function", use the function "methodFunction" as a way of combining genes
#    Alternatively, use one of the built-in functions 
#    Note: methodFunction must be a function that takes a vector of numbers as input and
#     outputs a single number. This function will return(0) or crash otherwise.

   recMethods = c("function","MaxMean","maxVariance","MinMean","absMinMean","absMaxMean");
   imethod = pmatch(method, recMethods);
        
   if (is.na(imethod)) 
      stop("Error: entered method is not a legal option. Recognized options are\n",
           "       *maxVariance*, *MaxMean*, *MinMean*, *absMaxMean*, *absMinMean*\n",
           "       or *function* for a user-defined function.")

   if (imethod > 1) 
   {
     selectionStatisticFnc = spaste(".cr.", method);
     selStatFnc = get(selectionStatisticFnc, mode = "function")
   } else {
     selStatFnc = match.fun(selectionStatisticFnc);
     if((!is.function(selStatFnc))&(!is.null(selStatFnc)))
            stop("Error: 'selectionStatisticFnc must be a function... please read the help file.")
   }
      
## Format the variables for use by this function
   colID[is.na(colID)] = group[is.na(colID)]    # Use group if row is missing
   rnDat[is.na(rnDat)]   = group[is.na(rnDat)];
   mdx = mtd.setColnames(mdx, rnDat);

   remove       = (is.na(colID))|(is.na(group)) # Omit if both gene and probe are missing
   colID  = colID[!remove];
   group = group[!remove];
   names(group) = colID
   colID = sort(intersect(rnDat,colID))
   if (length(colID)<=1)
      stop("None of the variable names in 'mdx' are in 'colID'.")

   group = group[colID]
   mdx  = mtd.apply(mdx, as.matrix);
   keepVars[keepVars] =  mtd.colnames(mdx) %in% colID;
   mdx = mtd.subset(mdx, , colID);

   probes = mtd.colnames(mdx)
   genes  = group[probes]
   tGenes = table(genes)
   colnames.out = sort(names(tGenes));
    
   if (getRepresentativeData)
   {
     mdxOut = mtd.apply(mdx, function(x) 
     {
       out = matrix(0, nrow(x), length(tGenes));
       rownames(out) = rownames(x);
       colnames(out) = colnames.out
       out;
     });
     names(mdxOut) = names(mdx);
   }

   representatives = rep("", length(colnames.out))
   names(representatives) = colnames.out; 
   
##  If !is.null(connectivityPower), default to the connectivity method with power=method
#      Collapse genes with multiple probe sets together using the following algorthim:
#      1) If there is one ps/g = keep
#      2) If there are 2 ps/g = (use "method" or "methodFunction")
#      3) If there are 3+ ps/g = take the max connectivity
#   Otherwise, use "method" if there are 3+ ps/g as well. 
   if(!is.null(connectivityPower)){
     if(!is.numeric(connectivityPower))
        stop("Error: if entered, connectivityPower must be numeric.")
     if(connectivityPower<=0)
       stop("Warning: connectivityPower must be >= 0.");

     if(any(nSamples<=5)){
       write("Warning: 5 or fewer samples, this method of probe collapse is unreliable...","")
       write("...Running anyway, but we suggest trying another method (for example, *mean*).","")
     }
   }
   
   # Run selectionStatisticFnc on all data; if quantile normalization is requested, normalize the selection
   # statistics across data sets.

   selectionStatistics = mtd.apply(mdx, function(x) 
             do.call(selStatFnc, c(list(x), statisticFncArguments)),
              mdaSimplify = TRUE);

   #if (FALSE) xxx = selectionStatistics;

   if (is.null(dim(selectionStatistics)))
      stop("Calculation of selection statistics produced results of zero or unqual lengths.");

   if (calibration=="full quantile")
      selectionStatistics = normalize.quantiles(selectionStatistics);

   #if (FALSE)
   #{
   #   sizeGrWindow(14, 5);
   #   par(mfrow = c(1,3));
   #   for (set in 1:nSets)
   #     hist(xxx[, set], breaks = 200);
#
##      for (set in 1:nSets)
 #       verboseScatterplot(xxx[, set], selectionStatistics[, set], samples = 10000)
 #  }

   consensusSelStat = .consensusCalculation(selectionStatistics, useMean = FALSE, setWeightMat = NULL,
                             consensusQuantile = consensusQuantile)$consensus;

   # Actually run the summarization.

   ones = sort(names(tGenes)[tGenes==1])
   if(useGroupHubs){
      twos = sort(names(tGenes)[tGenes==2]) # use "method" and connectivity
      more = sort(names(tGenes)[tGenes>2])
   } else { 
      twos = sort(names(tGenes)[tGenes>1]) # only use "method"
      more = character(0)
   }
   ones2genes =  match(ones, genes);
   if (getRepresentativeData) for (set in 1:nSets)
       mdxOut[[set]]$data[,ones] = mdx[[set]]$data[, ones2genes];
   representatives[ones] = probes[ones2genes];
   count = 0;

   if (length(twos) > 0)
   {
     if (verbose > 0)
       printFlush(spaste(spaces, "..selecting representatives for 2-variable groups.."));
     if (verbose > 1) pind = initProgInd(paste(spaces, ".."));
     repres = rep(NA, length(twos));
     for (ig in 1:length(twos))
     {
        g = twos[ig];
        probeIndex = which(genes==g);
        repres[ig] = probeIndex[which.max(consensusSelStat[probeIndex])];
        if (verbose > 1) pind = updateProgInd(ig/length(twos), pind);
     }
     if (verbose > 1) printFlush("");
     if (getRepresentativeData) for (set in 1:nSets)
        mdxOut[[set]]$data[, twos] = mdx[[set]]$data[, repres];
     representatives[twos] = probes[repres];
   }
   if (length(more) > 0)
   {
     if (verbose > 0)
       printFlush(spaste(spaces, "..selecting representatives for 3-variable groups.."));
     if (verbose > 1) pind = initProgInd(paste(spaces, ".."));
     genes.more = genes[genes %in% more];
     nAll = length(genes.more);
     connectivities = matrix(NA, nAll, nSets);
     for (ig in 1:length(more))
     {
        g = more[ig];
        keepProbes1 = which(genes==g);
        keep.inMore = which(genes.more==g);
        mdxTmp = mtd.subset(mdx, , keepProbes1);
        adj = mtd.apply(mdxTmp, function(x) do.call(adjacency, 
                c(list(x, type = "signed", power = connectivityPower), adjacencyArguments)));
        connectivities[keep.inMore, ] = mtd.apply(adj, colSums, mdaSimplify = TRUE);
        count = count + 1;
        if (count %% 50000 == 0) collectGarbage();
        if (verbose > 1) pind = updateProgInd(ig/(2*length(more)), pind);
     }

     if (calibration=="full quantile")
       connectivities = normalize.quantiles(connectivities);

     consConn = .consensusCalculation(connectivities, useMean = FALSE, setWeightMat = NULL,
                                               consensusQuantile = consensusQuantile)$consensus;
     repres.inMore = rep(0, length(more));
     for (ig in 1:length(more))
     {
        probeIndex = which(genes.more==more[ig]);
        repres.inMore[ig] = probeIndex[which.max(consConn[probeIndex])];
        if (verbose > 1) pind = updateProgInd(ig/(2*length(more)) + 0.5, pind);
     }
     repres = which(genes %in% more)[repres.inMore];
     if (verbose > 1) printFlush("");
     if (getRepresentativeData) for (set in 1:nSets)
       mdxOut[[set]]$data[, more] = as.numeric(mdx[[set]]$data[, repres]); 
     representatives[more] = probes[repres];
   }
      
   # Retreive the information about which probes were saved, and include that information
   #   as part of the output.  

   out2 = cbind(colnames.out, representatives)
   colnames(out2) = c("group","selectedColID")

   reprIndicator = keepVars;
   reprIndicator[keepVars] [match(representatives, mtd.colnames(mdx))] = TRUE
   reprIndicator = colnames.in
   out = list(representatives = out2, 
              varSelected = reprIndicator,
              representativeData = if (getRepresentativeData) mdxOut else NULL)
   return(out)
      
} 

