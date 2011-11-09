
# The function standardScreeningBinaryTrait computes widely used statistics for relating the columns of
# the input data frame (argument datExpr) to a binary sample trait (argument y). The statistics include
# Student t-test p-value and the corresponding local false discovery rate (known as q-value, Storey et al
# 2004), the fold change, the area under the ROC curve (also known as C-index), mean values etc. If the
# input option kruskalTest is set to TRUE, it also computes the kruskal Wallist test p-value and
# corresponding q-value. The kruskal Wallis test is a non-parametric, rank-based group comparison test.

standardScreeningBinaryTrait=function(datExpr, y, 
           corFnc = cor, corOptions = list(use = 'p'),
           kruskalTest=FALSE, qValues = FALSE, var.equal=FALSE, na.action="na.exclude") 
{
   datExpr=data.frame(datExpr)
   levelsy=levels(factor(y))
   if (length(levelsy)>2 ) 
     stop("The sample trait y contains more than 2 levels. Please input a binary variable y")
   if (length(levelsy)==1 ) 
     stop("The sample trait y is constant. Please input a binary sample trait with some variation.")
   yNumeric=ifelse( y==levelsy[[1]], 1, ifelse( y==levelsy[[2]], 2, NA))
   if (length(yNumeric) !=dim(datExpr)[[1]] ) 
     stop("the length of the sample trait y does not equal the number of rows of datExpr")
   corPearson=rep(NA, dim(datExpr)[[2]] ) 
   pvalueStudent = t.Student = Z.Student = rep(NA, dim(datExpr)[[2]] ) 
   pvaluekruskal = stat.Kruskal = Z.Kruskal = sign.Kruskal = rep(NA, dim(datExpr)[[2]] ) 
   nPresent = rep(0, dim(datExpr)[[2]] ) 
   AreaUnderROC=rep(NA, dim(datExpr)[[2]] ) 

  if (var.equal) 
     printFlush(paste("Warning: T-test that assumes equal variances in each group is requested.\n",
                      "This is not the default option for t.test. We recommend to use var.equal=FALSE."));

  corFnc = match.fun(corFnc);
  corOptions$y = yNumeric;
  for (i in 1:dim(datExpr)[[2]]) {
        corOptions$x = as.numeric(datExpr[,i]);
        corPearson[i] = as.numeric(do.call(corFnc, corOptions));
        no.present=  sum( ! is.na(datExpr[,i])  & ! is.na(yNumeric)   )
        nPresent[i] = no.present;
        no.present1=  sum( ! is.na(datExpr[,i])  & !is.na(yNumeric) & yNumeric==1  )
        no.present2=  sum( ! is.na(datExpr[,i])  & !is.na(yNumeric) & yNumeric==2  )
        if (no.present1<2 | no.present2<2 ) 
        {
           pvalueStudent[i]= t.Student[i] = NA 
        } else {
          tst = t.test( as.numeric(datExpr[,i])~yNumeric,var.equal=var.equal,na.action=na.action)
          pvalueStudent[i] = tst$p.value;
          t.Student[i] = -tst$statistic 
          # The - sign above is intentional to make the sign of t consistent with correlation
        }
        AreaUnderROC[i] = rcorr.cens(datExpr[, i], yNumeric, outx = TRUE)[[1]]
        if (kruskalTest) {
            if (no.present<5 ) 
            {
               pvaluekruskal[i] = stat.Kruskal[i] = NA
            } else {
               kt = kruskal.test(datExpr[, i] ~ factor(yNumeric),  na.action="na.exclude")
               pvaluekruskal[i] = kt$p.value;
               stat.Kruskal[i] = kt$statistic;
               # Find which side is higher
               r = rank(datExpr[, i]);
               sums = tapply(r, factor(yNumeric), sum);
               sign.Kruskal[i] = 2 * ( (sums[1] < sums[2]) - 0.5);
               # sign.Kruskal is 1 if the ranks in group 1 are smaller than in group 2
           }
        } 
    }
    q.Student=rep(NA, length(pvalueStudent) )
    rest1= ! is.na(pvalueStudent) 
    if (qValues) 
    {
      x = try({ q.Student[rest1] = qvalue(pvalueStudent[rest1])$qvalues }, silent = TRUE);
      if (inherits(x, "try-error"))
        printFlush(paste("Warning in standardScreeningBinaryTrait: function qvalue returned an error.\n",
                         "calculated q-values will be invalid. qvalue error:\n\n", x, "\n")) 
      if (kruskalTest) {
         q.kruskal=rep(NA, length(pvaluekruskal) )
         rest1= ! is.na(pvaluekruskal) 
         xx = try( { q.kruskal[rest1] = qvalue(pvaluekruskal[rest1])$qvalues} , silent = TRUE);
         if (inherits(xx, "try-error"))
           printFlush(paste("Warning in standardScreeningBinaryTrait: function qvalue returned an error.\n",
                            "calculated q-values will be invalid. qvalue error:\n\n", xx, "\n")) 
      }
    }
    meanLevel1 = as.numeric(apply(datExpr[y == levelsy[[1]] & !is.na(y), ], 2, mean, na.rm = TRUE));
    meanLevel2 = as.numeric(apply(datExpr[y == levelsy[[2]] & !is.na(y), ], 2, mean, na.rm = TRUE));
  
    Z.Student = qnorm(pvalueStudent/2, lower.tail = FALSE) * sign(t.Student);
    if (kruskalTest)
      Z.Kruskal = qnorm(pvaluekruskal/2, lower.tail = FALSE) * sign(stat.Kruskal);

    
    stderr1=function(x) {no.present=sum(!is.na(x));
    if (no.present<2) out=NA else {out=sqrt(var(x,na.rm=TRUE)/no.present) }
    out } # end of function stderr1

    SE.Level1 = as.numeric(apply(datExpr[y == levelsy[[1]] & !is.na(y), ], 2, stderr1))
    SE.Level2 =as.numeric(apply(datExpr[y == levelsy[[2]] & !is.na(y), ], 2, stderr1))

    FoldChangeLevel1vsLevel2 = ifelse(meanLevel1/meanLevel2 > 
        1, meanLevel1/meanLevel2, -meanLevel2/meanLevel1)
    
    output = data.frame(ID = dimnames(datExpr)[[2]], corPearson = corPearson, 
        t.Student = t.Student,
        pvalueStudent = pvalueStudent, 
        FoldChange = FoldChangeLevel1vsLevel2, 
         meanFirstGroup = meanLevel1,
        meanSecondGroup = meanLevel2, 
         SE.FirstGroup = SE.Level1,
      SE.SecondGroup = SE.Level2,
      AreaUnderROC = AreaUnderROC)

    if (kruskalTest) {
        output = data.frame(output, stat.Kruskal = stat.Kruskal, 
                            stat.Kruskal.signed = sign.Kruskal * stat.Kruskal,
                            pvaluekruskal = pvaluekruskal);
    }

   if (qValues && !inherits(x, "try-error")) output=data.frame(output, q.Student)
   if (qValues &  kruskalTest ) {
      if ( !inherits(xx, "try-error")) output=data.frame(output, q.kruskal)
   }
   names(output)[3:5] = paste(names(output)[3:5], levelsy[[1]], "vs", levelsy[[2]], sep = ".")
   output = data.frame(output, nPresentSamples = nPresent);
   output
}




