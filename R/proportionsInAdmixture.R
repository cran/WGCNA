proportionsInAdmixture<-function 
(MarkerMeansPure, 
datE.Admixture, 
calculateConditionNumber = FALSE,
coefToProportion = TRUE)
{
    datE.Admixture = data.frame(datE.Admixture)
    if (sum(is.na(names(datE.Admixture))) > 0) {
        warning("Some of the column names of datE.Admixture are missing. Recommendation: check or assign column names. But for your convenience, we remove the corresponding columns of datE.Admixture from the analysis.")
        datE.Admixture = datE.Admixture[, !is.na(names(datE.Admixture))]
    }
    if (sum(names(datE.Admixture) == "") > 0) {
        warning("Some of the column names of datE.Admixture are missing. Recommendation: check or assign column names.  But for your convenience, we remove the corresponding columns of datE.Admixture from the analysis.")
        datE.Admixture = datE.Admixture[, names(datE.Admixture) != 
            ""]
    }
    MarkerID = MarkerMeansPure[, 1]
    if (sum(is.na(MarkerID)) > 0) {
        warning("Some of the marker are missing (NA). Recommendation: check the first column of the input MarkerMeansPure. It should contain marker names. But for your convenience, we remove the corresponding markers from the analysis.")
        MarkerMeansPure = MarkerMeansPure[!is.na(MarkerID), ]
        MarkerID = MarkerMeansPure[, 1]
    }
    if (sum(MarkerID == "", na.rm = T) > 0) {
        warning("Some of the marker names are empty strings. Recommendation: check the first column of the input MarkerMeansPure. It should contain marker names.  But for your convenience, we remove the corresponding markers from the analysis.")
        MarkerMeansPure = MarkerMeansPure[MarkerID != "", ]
        MarkerID = MarkerMeansPure[, 1]
    }
    noMissingValuesMarker = as.numeric(apply(is.na(MarkerMeansPure[, 
        -1]), 1, sum))
    if (max(noMissingValuesMarker, na.rm = T) > 0) {
        warning("Some of the markers  (rows of MarkerMeansPure) contain missing values. This is problematic.\nFor your convenience, we remove the corresponding markers (rows) from the analysis.")
        MarkerMeansPure = MarkerMeansPure[noMissingValuesMarker == 
            0, ]
        MarkerID = MarkerMeansPure[, 1]
    }
    match1 = match(MarkerID, names(datE.Admixture))
    match1 = match1[!is.na(match1)]
    if (length(match1) == 0) 
        stop("None of the marker names correspond to column names of the input datE.Admixture. Possible solutions: Transpose datE.Admixture or MarkerMeansPure. Or make sure to assign suitable names to the columns of datE.Admixture, e.g. as follows dimnames(datE.Admixture)[[2]]=GeneSymbols.")
    if (length(match1) < dim(MarkerMeansPure)[[1]]) {
        warning(paste("Only", length(match1), "out of ", dim(MarkerMeansPure)[[1]], 
            "rows of MarkerMeansPure correspond to columns of datE.Admixture. \nIf this suprises you, check the the first column of MarkerMeansPure or the column names of datE.Admixture. \nThe output contains a list of markers that could be identified."))
    }
    datE.MarkersAdmixtureTranspose = t(datE.Admixture[, match1])
    match2 = match(names(datE.Admixture)[match1], MarkerID)
    match2 = match2[!is.na(match2)]
    MarkerMeansPure = MarkerMeansPure[match2, ]
    if (sum(as.character(MarkerMeansPure[, 1]) != dimnames(datE.MarkersAdmixtureTranspose)[[1]], 
        na.rm = T) > 0) 
        stop("I am sorry but things do not line up. Maybe you need to look inside the R code. Specifically,\nas.character(MarkerMeansPure) != dimnames(datE.MarkersAdmixtureTranspose)[[1]]")
    conditionNumber = NA
    if (dim(MarkerMeansPure)[[2]] == 2) {
        A = as.matrix(MarkerMeansPure[, -1], ncol = 1)
    }
    else {
        A = as.matrix(MarkerMeansPure[, -1])
    }
    if (dim(as.matrix(A))[[2]] > 1 & dim(as.matrix(A))[[1]] > 
        1 & calculateConditionNumber) {
        conditionNumber = kappa(A)
    }
    datCoef = t(lm(datE.MarkersAdmixtureTranspose ~ A)$coefficients[-1, 
        ])
    coef2prop = function(coef) {
        prop = rep(NA, length(coef))
        coef[coef < 0] = 0
        if (sum(coef, na.rm = T) > 0 & !is.na(sum(coef, na.rm = T))) {
            prop = coef/sum(coef, na.rm = T)
        }
        prop
    }
    if (coefToProportion) {
        PredictedProportions = data.frame(t(apply(datCoef, 1, 
            coef2prop)))
    }
    else {
        PredictedProportions = datCoef
    }
    dimnames(PredictedProportions)[[1]] = dimnames(datE.Admixture)[[1]]
    out = list(PredictedProportions = PredictedProportions, datCoef = datCoef, 
        conditionNumber = conditionNumber, markersUsed = as.character(MarkerMeansPure[, 
            1]))
    out
}
