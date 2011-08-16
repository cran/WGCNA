#
userListEnrichment <- function (geneR, labelR, fnIn = NULL, catNmIn = fnIn, nameOut = "enrichment.csv", 
    useBrainLists = FALSE, useBloodAtlases=FALSE, omitCategories = "grey", outputCorrectedPvalues = TRUE) 
{
    if (length(geneR) != length(labelR)) 
        stop("geneR and labelR must have same number of elements.")
    if (length(catNmIn) < length(fnIn)) {
        catNmIn = c(catNmIn, fnIn[(length(catNmIn) + 1):length(fnIn)])
        write("WARNING: not enough category names.  \n\t\t\t   Naming remaining categories with file names.", 
            "")
    }
    if (is.null(fnIn) & (! (useBrainLists | useBloodAtlases)) ) 
        stop("Either enter user-defined lists or set useBrainLists or useBloodAtlases = TRUE.")

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
    if(useBloodAtlases) { 
        if (!(exists("BloodLists"))) BloodLists = NULL;
	 	data("BloodLists",envir=sys.frame(sys.nframe()));
		write("See help file for details regarding blood atlas references.",
		"")
	glIn = rbind(glIn, cbind(BloodLists, Type = rep("Blood", nrow(BloodLists))))
    }
    removeDups = unique(paste(as.character(glIn[, 1]), as.character(glIn[, 
        2]), sep = "@#$%"))
    if (length(removeDups) < length(glIn[, 1])) 
        glIn = t(as.matrix(as.data.frame(strsplit(removeDups, "@#$%", fixed = TRUE))))
    geneIn = as.character(glIn[, 1])
    labelIn = as.character(glIn[, 2])
    geneAll = sort(unique(geneR))
    keep = is.element(geneIn, geneAll)
    geneIn = geneIn[keep]
    labelIn = labelIn[keep]
    catsR = sort(unique(labelR))
    catsR = catsR[!is.element(catsR, omitCategories)]
    catsIn = sort(unique(labelIn))
    typeIn = glIn[keep, ][match(catsIn, labelIn), 3];
    lenAll = length(geneAll)
    results = list(pValues = NULL, ovGenes = list(), sigOverlaps = NULL)
    namesOv = NULL
    for (r in 1:length(catsR)) for (i in 1:length(catsIn)) {
        gr = geneR[labelR == catsR[r]]
        gi = geneIn[labelIn == catsIn[i]]
        go = intersect(gr, gi)
        lr = length(gr)
        li = length(gi)
        lo = length(go)
        pv = .phyper2(lenAll, lr, li, lo, FALSE)
        if (lo == 0) 
            pv = 1
        if (pv < 1e-04) 
            pv = .phyper2(lenAll, lr, li, lo, TRUE)
        pOut = c(catsR[r], catsIn[i], typeIn[i], lo, pv)
        results$pValues = rbind(results$pValues, pOut)
        namesOv = c(namesOv, paste(catsR[r], "--", catsIn[i]))
        results$ovGenes[[length(namesOv)]] = go
    }
    results$pValues = cbind(results$pValues, apply(cbind(1, as.numeric(results$pValues[, 
        4]) * length(namesOv)), 1, min))
    colnames(results$pValues) = c("InputCategories", "UserDefinedCategories", "Type", 
        "NumOverlap", "Pvalues", "CorrectedPvalues")
    names(results$ovGenes) = namesOv
    results$sigOverlaps = results$pValues[as.numeric(results$pValues[, 
        5]) < 0.05, c(1, 2, 3, 5)]
    if (!outputCorrectedPvalues) {
        results$sigOverlaps = results$pValues[as.numeric(results$pValues[, 
            4]) < 0.05, c(1, 2, 3, 4)]
        write("Note that outputted p-values are not corrected for multiple comparisons.", 
            "")
    }
    results$sigOverlaps = as.data.frame(results$sigOverlaps[order(as.numeric(results$sigOverlaps[, 
        4])), ]);
    row.names(results$sigOverlaps) = NULL;
    write.csv(results$sigOverlaps, nameOut, row.names = FALSE)
    write(paste(length(namesOv), "comparisons were successfully performed."), 
        "")
    return(results)
}
