# Accuracy measures, modified from the WGCNA version.

# Helper function: contingency table of 2 variables that will also include rows/columns for levels that do
# not appear in x or y.

.table2.allLevels = function(x, y, levels.x = sort(unique(x)), levels.y = sort(unique(y)), setNames = FALSE)
{
  nx = length(levels.x);
  ny = length(levels.y);
  t = table(x, y);

  out = matrix(0, nx, ny);
  if (setNames)
  {
    rownames(out) = levels.x;
    colnames(out) = levels.y;
  }
  out[ match(rownames(t), levels.x), match(colnames(t), levels.y) ] = t;
  out;
}

# accuracy measures

accuracyMeasures = function(predicted, observed = NULL, 
          type = c("auto", "binary", "quantitative"),
          levels = if (isTRUE(all.equal(dim(predicted), c(2,2)))) colnames(predicted)
                     else if (is.factor(predicted))
                        sort(unique(c(as.character(predicted), as.character(observed))))
                     else sort(unique(c(observed, predicted))),
          negativeLevel = levels[2], positiveLevel = levels[1] )
{
  type = match.arg(type);
  if (type=="auto")
  {
    if (!is.null(dim(predicted)))
    {
      if (isTRUE(all.equal(dim(predicted), c(2,2))))
      {
        type = "binary"
      } else
        stop("If supplying a matrix in 'predicted', it must be a 2x2 contingency table.");
    } else {
      if (is.null(observed)) 
        stop("When 'predicted' is a vector, 'observed' must be given and have the same length as 'predicted'.");

      if (length(levels)==2) 
      {
        type = "binary"
      } else
        type = "quantitative"
    }
  }

  if (type=="binary")
  {
    if (is.null(dim(predicted)))
    {
      if (is.null(observed)) 
        stop("When 'predicted' is a vector, 'observed' must be given and have the same length as 'predicted'.");
      if ( length(predicted)!=length(observed) )
        stop("When both 'predicted' and 'observed' are given, they must be vectors of the same length.");
      if (length(levels)!=2) 
        stop("'levels' must contain 2 entries (the possible values of the binary variables\n", 
             "   'predicted' and 'observed').");

      tab = .table2.allLevels(predicted, observed, levels.x = levels, levels.y = levels, setNames = TRUE);
    } else {
      tab = predicted;
      if (is.null(colnames(tab)) | is.null(rownames(tab)))
        stop("When 'predicted' is a contingency table, it must have valid colnames and rownames.");

    }

    if (  ncol(tab) !=2 |  nrow(tab) !=2 ) 
      stop("The input table must be a 2x2 table. ")

    if (negativeLevel==positiveLevel) 
      stop("'negativeLevel' and 'positiveLevel' cannot be the same.");

    neg = match(negativeLevel, colnames(tab));
    if (is.na(neg))
      stop(spaste("Cannot find the negative level ", negativeLevel, 
                  " among the colnames of the contingency table.\n   Please check the input and try again."))
    pos = match(positiveLevel, colnames(tab));
    if (is.na(pos))
      stop(spaste("Cannot find the positive level ", positiveLevel, 
                  " among the colnames of the contingency table.\n   Please check the input and try again."))
      
    if (  sum(is.na(tab) ) ) 
      warning("Missing data should not be present in input.\n", 
              "  Suggestion: check whether NA should be coded as 0.")

    is.wholenumber =function(x, tol = .Machine$double.eps^0.5) { abs(x - round(x)) < tol }

    if (  sum( !is.wholenumber(tab), na.rm=T  ) >0) 
      warning("STRONG WARNING: The input table contains non-integers, which does not make sense.")

    if (  sum( tab<0, na.rm=T  ) >0) 
      stop("The input table cannot contain negative numbers.");

    num1=sum(diag(tab),na.rm=T)
    denom1=sum(tab,na.rm=T)
    if (denom1==0) 
       warning("The input table has zero observations (sum of all cells is zero).")

    TP=tab[pos, pos]
    FP=tab[pos, neg]
    FN=tab[neg, pos]
    TN=tab[neg, neg]

    error.rate= ifelse(denom1==0,NA, 1-num1/denom1)
    Accuracy= ifelse(denom1==0,NA,  num1/denom1 )
    Specificity= ifelse(FP + TN==0, NA,  TN / (FP + TN) )
    Sensitivity= ifelse(TP + FN==0, NA,  TP / (TP + FN) )
    NegativePredictiveValue= ifelse(FN + TN==0,NA,  TN / (FN + TN) )
    PositivePredictiveValue=ifelse(TP + FP==0,NA,    TP / (TP + FP) )
    FalsePositiveRate = 1 - Specificity 
    FalseNegativeRate = 1 - Sensitivity 
    Power = Sensitivity 
    LikelihoodRatioPositive = ifelse(1 - Specificity==0,NA, Sensitivity / (1 - Specificity) )
    LikelihoodRatioNegative = ifelse(Specificity==0, NA,  (1 - Sensitivity) / Specificity )
    NaiveErrorRate = ifelse(denom1==0,NA,   
                             min(c(tab[pos, pos]+ tab[neg, pos] , tab[pos, neg]+ tab[neg, neg] ))/denom1   ) 
    out=data.frame(
          Measure= c("Error.Rate","Accuracy", "Specificity","Sensitivity","NegativePredictiveValue",
                     "PositivePredictiveValue","FalsePositiveRate","FalseNegativeRate","Power",
                     "LikelihoodRatioPositive","LikelihoodRatioNegative", "NaiveErrorRate", "NegativeLevel",
                     "PositiveLevel"),
          Value=c(error.rate,Accuracy, Specificity,Sensitivity,NegativePredictiveValue,
                  PositivePredictiveValue,FalsePositiveRate,FalseNegativeRate,Power,
                  LikelihoodRatioPositive,LikelihoodRatioNegative,NaiveErrorRate, negativeLevel,
                  positiveLevel));
  } else if (type=="quantitative") 
  {
     if (!is.null(dim(predicted))) 
       stop("When 'type' is \"quantitative\", 'predicted' cannot be a 2-dimensional matrix.");
     if (length(predicted)!=length(observed))
       stop("'predicted' and 'observed' must be vectors of the same length.");

     cr = cor(predicted, observed, use = 'p');
     out = data.frame(
          Measure = c("Cor", "R.squared", "MeanSquareError", "MedianAbsoluteError", "Cindex"),
          Value = c(cr, 
                    cr^2, 
                    mean( (predicted-observed)^2,na.rm=TRUE), 
                    median((predicted-observed)^2,na.rm=TRUE),
                    rcorr.cens(predicted,observed,outx=TRUE)[[1]]));
  }

  out;
}
