#
userListEnrichment <- function (geneR, labelR, fnIn = NULL, catNmIn = fnIn, nameOut = "enrichment.csv", 
    useBrainLists = FALSE, useBloodAtlases = FALSE, omitCategories = "grey", 
    outputCorrectedPvalues = TRUE, useStemCellLists = FALSE, outputGenes = FALSE, 
  	minGenesInCategory = 1, useBrainRegionMarkers = FALSE, useImmunePathwayLists = FALSE,
	usePalazzoloWang = FALSE) 
{
    if (length(geneR) != length(labelR)) 
        stop("geneR and labelR must have same number of elements.")
    if (length(catNmIn) < length(fnIn)) {
        catNmIn = c(catNmIn, fnIn[(length(catNmIn) + 1):length(fnIn)])
        write("WARNING: not enough category names.  \n\t\t\t   Naming remaining categories with file names.", 
            "")
    }
    if (is.null(fnIn) & (! (useBrainLists | useBloodAtlases | useStemCellLists | useBrainRegionMarkers
	                      | useImmunePathwayLists | usePalazzoloWang)) ) 
        stop("Either enter user-defined lists or set one of the use_____ parameters to TRUE.")

    glIn = NULL
    if (length(fnIn)>0) 
    {
      for (i in 1:length(fnIn)) 
      {
        ext = substr(fnIn[i], nchar(fnIn[i]) - 2, nchar(fnIn[i]))
        if (ext == "csv") {
            datIn = read.csv(fnIn[i])
            if (colnames(datIn)[2] == "Gene") {
                datIn = datIn[, 2:3]
            }
            else {
                datIn = datIn[, 1:2]
            }
        }
        else {
            datIn = scan(fnIn[i], what = "character", sep = "\n")
            datIn = cbind(datIn[2:length(datIn)], datIn[1])
        }
        colnames(datIn) = c("Gene", "Category")
        datIn[, 2] = paste(datIn[, 2], catNmIn[i], sep = "__")
        glIn = rbind(glIn, datIn)
      }
      glIn = cbind(glIn, Type = rep("User", nrow(glIn)));
    }
    if (useBrainLists) {
        if (!(exists("BrainLists"))) BrainLists = NULL;
        data("BrainLists",envir=sys.frame(sys.nframe()));
        write("See help file for details regarding brain list references.", 
            "")
        glIn = rbind(glIn, cbind(BrainLists, Type = rep("Brain", nrow(BrainLists))))
    }
    if (useBloodAtlases) { 
        if (!(exists("BloodLists"))) BloodLists = NULL;
        data("BloodLists",envir=sys.frame(sys.nframe()));
        write("See help file for details regarding blood atlas references.",
		        "")
	      glIn = rbind(glIn, cbind(BloodLists, Type = rep("Blood", nrow(BloodLists))))
    }
    if (useStemCellLists) {
        if (!(exists("SCsLists"))) SCsLists = NULL;
        data("SCsLists",envir=sys.frame(sys.nframe()));
        write("See help file for details regarding stem cell list references.", 
            "")
        glIn = rbind(glIn, cbind(SCsLists, Type = rep("StemCells", nrow(SCsLists))))
    }
	  if (useBrainRegionMarkers) {
        if (!(exists("BrainRegionMarkers"))) BrainRegionMarkers = NULL;
        data("BrainRegionMarkers",envir=sys.frame(sys.nframe()));
        write("Brain region markers from http://human.brain-map.org/ -- see help file for details.", 
            "")
        glIn = rbind(glIn, cbind(BrainRegionMarkers, Type = rep("HumanBrainRegions", nrow(BrainRegionMarkers))))
    }
    if (useImmunePathwayLists) { 
        if (!(exists("ImmunePathwayLists"))) ImmunePathwayLists = NULL;
        data("ImmunePathwayLists",envir=sys.frame(sys.nframe()));
        write("See help file for details regarding immune pathways.",
		        "")
	      glIn = rbind(glIn, cbind(ImmunePathwayLists, Type = rep("Immune", nrow(ImmunePathwayLists))))
    }
    if (usePalazzoloWang) { 
        if (!(exists("PWLists"))) PWLists = NULL;
        data("PWLists",envir=sys.frame(sys.nframe()));
        write("See help file for details regarding Palazzolo / Wang lists from CHDI.",
		        "")
		write("---- there are many of these gene sets so the function may take several minutes to run.",
		        "")
	      glIn = rbind(glIn, cbind(PWLists, Type = rep("PW_Lists", nrow(PWLists))))
    }

    #removeDups = unique(paste(as.character(glIn[, 1]), as.character(glIn[, 
    #    2]), as.character(glIn[, 3]), sep = "@#$%"))
    #if (length(removeDups) < length(glIn[, 1])) 
    #    glIn = t(as.matrix(as.data.frame(strsplit(removeDups, "@#$%", fixed = TRUE))))

    
    glIn.2 = glIn[!duplicated(as.data.frame(glIn)), ]

    geneIn = as.character(glIn.2[, 1])
    labelIn = as.character(glIn.2[, 2])
    geneAll = sort(unique(geneR))
    keep = is.element(geneIn, geneAll)
    geneIn = geneIn[keep]
    labelIn = labelIn[keep]
    catsR = sort(unique(labelR))
    omitCategories = c(omitCategories, "background")      
    catsR = catsR[!is.element(catsR, omitCategories)]
    catsIn = sort(unique(labelIn))
    typeIn = glIn.2[keep, ][match(catsIn, labelIn), 3];
    lenAll = length(geneAll)
    nCols.pValues = 5;
    nComparisons = length(catsR) * length(catsIn);
    nIn = length(catsIn);
    nR = length(catsR);
    index = 1; 
    nOverlap = rep(0, nComparisons);
    pValues = rep(1, nComparisons);
    ovGenes = vector(mode = "list", length = nComparisons);
    isI = matrix(FALSE, lenAll, nIn);
    for (i in 1:nIn)
    {
      isI[, i] = is.element(geneAll,geneIn[labelIn == catsIn[i]]);
    }
    for (r in 1:length(catsR)) 
    {
      isR  = is.element(geneAll,geneR[(labelR == catsR[r])])
      for (i in 1:length(catsIn)) 
      { 
        isI.1 = isI[, i];
        lyn  = sum(isR&(!isI.1))
        lny  = sum(isI.1&(!isR))
        lyy  = sum(isR&isI.1)
        gyy  = geneAll[isR&isI.1]
        lnn  = lenAll - lyy - lyn - lny
        pv   = fisher.test(matrix(c(lnn,lny,lyn,lyy), 2, 2), alternative = "greater")$p.value
        nOverlap[index] = lyy;
        pValues[index] = pv;
        ovGenes[[index]] = gyy
        index = index + 1
      }
    }

    results = list(pValues = data.frame(InputCategories = rep(catsR, rep(nIn, nR)),
                                        UserDefinedCategories = rep(catsIn, nR),
                                        Type = rep(typeIn, nR),
                                        NumOverlap = nOverlap,
                                        Pvalues = pValues,
                                        CorrectedPvalues = ifelse(pValues * nComparisons > 1, 1,
                                                                  pValues * nComparisons)),
                   ovGenes = ovGenes);
    namesOv = paste(results$pValues$InputCategories, "--", results$pValues$UserDefinedCategories);
    names(results$ovGenes) = namesOv
    if (outputCorrectedPvalues) {
        results$sigOverlaps = results$pValues[results$pValues$CorrectedPvalues < 0.05, c(1, 2, 3, 6)]
    } else {
        results$sigOverlaps = results$pValues[results$pValues$Pvalues < 0.05, c(1, 2, 3, 5)]
        write("Note that outputted p-values are not corrected for multiple comparisons.", 
            "")
    }
    results$sigOverlaps = results$sigOverlaps[order(results$sigOverlaps[, 4]), ];
    row.names(results$sigOverlaps) = NULL;
    
    rSig  = results$sigOverlaps
    nSig = nrow(rSig);
    if (nSig > 0)
    {
       rCats = paste(rSig$InputCategories,"--",rSig$UserDefinedCategories)
       rNums <- rep(0, nSig);
       rGenes <- rep("", nSig);
       for (i in 1:nSig)
       {
         rGn    = results$ovGenes[[which(names(results$ovGenes)==rCats[i])]]
         rNums[i]  = length(rGn);
         rGenes[i] = paste(rGn,collapse=", ");
       }
       rSig$NumGenes = rNums
       rSig$CategoryGenes = rGenes
       rSig = rSig[rSig$NumGenes>=minGenesInCategory,]
       if(!outputGenes) rSig = rSig[,1:4]
       results$sigOverlaps = rSig
       write.csv(results$sigOverlaps, file = nameOut, row.names = FALSE)
    }

    write(paste(length(namesOv), "comparisons were successfully performed."), 
        "")
    return(results)
}
