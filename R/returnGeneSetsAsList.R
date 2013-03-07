returnGeneSetsAsList <- function (fnIn = NULL, catNmIn = fnIn, useBrainLists = FALSE, useBloodAtlases = FALSE, 
    useStemCellLists = FALSE, useBrainRegionMarkers = FALSE, useImmunePathwayLists = FALSE, geneSubset=NULL) 
{
    if (length(catNmIn) < length(fnIn)) {
        catNmIn = c(catNmIn, fnIn[(length(catNmIn) + 1):length(fnIn)])
        write("WARNING: not enough category names.  \n\t\t\t   Naming remaining categories with file names.", 
            "")
    }
    if (is.null(fnIn) & (!(useBrainLists | useBloodAtlases | 
        useStemCellLists | useBrainRegionMarkers | useImmunePathwayLists))) 
        stop("Either enter user-defined lists or set one of the use_____ parameters to TRUE.")
    glIn = NULL
    if (length(fnIn) > 0) {
        for (i in 1:length(fnIn)) {
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
        glIn = cbind(glIn, Type = rep("User", nrow(glIn)))
    }
    if (useBrainLists) {
        if (!(exists("BrainLists"))) 
            BrainLists = NULL
        data("BrainLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding brain list references.", 
            "")
        glIn = rbind(glIn, cbind(BrainLists, Type = rep("Brain", 
            nrow(BrainLists))))
    }
    if (useBloodAtlases) {
        if (!(exists("BloodLists"))) 
            BloodLists = NULL
        data("BloodLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding blood atlas references.", 
            "")
        glIn = rbind(glIn, cbind(BloodLists, Type = rep("Blood", 
            nrow(BloodLists))))
    }
    if (useStemCellLists) {
        if (!(exists("SCsLists"))) 
            SCsLists = NULL
        data("SCsLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding stem cell list references.", 
            "")
        glIn = rbind(glIn, cbind(SCsLists, Type = rep("StemCells", 
            nrow(SCsLists))))
    }
    if (useBrainRegionMarkers) {
        if (!(exists("BrainRegionMarkers"))) 
            BrainRegionMarkers = NULL
        data("BrainRegionMarkers", envir = sys.frame(sys.nframe()))
        write("Brain region markers from http://human.brain-map.org/ -- See userListEnrichment help file for details.", 
            "")
        glIn = rbind(glIn, cbind(BrainRegionMarkers, Type = rep("HumanBrainRegions", 
            nrow(BrainRegionMarkers))))
    }
    if (useImmunePathwayLists) {
        if (!(exists("ImmunePathwayLists"))) 
            ImmunePathwayLists = NULL
        data("ImmunePathwayLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding immune pathways.", 
            "")
        glIn = rbind(glIn, cbind(ImmunePathwayLists, Type = rep("Immune", 
            nrow(ImmunePathwayLists))))
    }
    removeDups = unique(paste(as.character(glIn[, 1]), as.character(glIn[, 
        2]), as.character(glIn[, 3]), sep = "@#$%"))
    if (length(removeDups) < length(glIn[, 1])) 
        glIn = t(as.matrix(as.data.frame(strsplit(removeDups, 
            "@#$%", fixed = TRUE))))
    geneIn = as.character(glIn[, 1])
    labelIn = paste(as.character(glIn[, 2]),as.character(glIn[, 3]),sep="__")
	if(!is.null(geneSubset)){
      keep = is.element(geneIn, geneSubset)
      geneIn = geneIn[keep]
      labelIn = labelIn[keep]
	}
	if(length(geneIn)<2)
	   stop("Please include a larger geneSubset, or set geneSubset=NULL.")
    allLabels <- sort(unique(labelIn))
	geneSet <- list()
	for (i in 1:length(allLabels))  geneSet[[i]] = geneIn[labelIn==allLabels[i]]
	names(geneSet) = allLabels
    return(geneSet)
}