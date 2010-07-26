# The function standardScreeningBinaryTrait computes widely used statistics for relating the columns of
# the input data frame (argument datExpr) to a binary sample trait (argument y). The statistics include
# Student t-test p-value and the corresponding local false discovery rate (known as q-value, Storey et al
# 2004), the fold change, the area under the ROC curve (also known as C-index), mean values etc. If the
# input option kruskalTest is set to TRUE, it also computes the kruskal Wallist test p-value and
# corresponding q-value. The kruskal Wallis test is a non-parametric, rank-based group comparison test.

standardScreeningBinaryTrait=function(datExpr, y, kruskalTest=FALSE) 
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
   pvalueStudent=rep(NA, dim(datExpr)[[2]] ) 
   pvaluekruskal=rep(NA, dim(datExpr)[[2]] ) 
   AreaUnderROC=rep(NA, dim(datExpr)[[2]] ) 

  for (i in 1:dim(datExpr)[[2]]) {
        corPearson[i] = as.numeric(cor(yNumeric, datExpr[,i], use = "p"))
        no.present=  sum( ! is.na(datExpr[,i])  & ! is.na(yNumeric)   )
        pvalueStudent[i] = corPvalueStudent(corPearson[i], no.present )
        AreaUnderROC[i] = rcorr.cens(datExpr[, i], yNumeric, outx = T)[[1]]
        if (kruskalTest) {
            if (no.present<5 )  pvaluekruskal[i]=NA
                   if (no.present>=5 )  {     
                      pvaluekruskal[i] = kruskal.test(datExpr[, i] ~ factor(yNumeric),  
                                                      na.action="na.exclude")$p.value
                   }
        } 
    }
    q.Student=rep(NA, length(pvalueStudent) )
    rest1= ! is.na(pvalueStudent) 
    q.Student[rest1] = qvalue(pvalueStudent[rest1])$qvalues
    if (kruskalTest) {
         q.kruskal=rep(NA, length(pvaluekruskal) )
                   rest1= ! is.na(pvaluekruskal) 
          q.kruskal[rest1] = qvalue(pvaluekruskal[rest1])$qvalues
    }
    meanLevel1 = apply(datExpr[y == levelsy[[1]] & !is.na(y), 
        ], 2, mean, na.rm = T)
    meanLevel2 = apply(datExpr[y == levelsy[[2]] & !is.na(y), 
        ], 2, mean, na.rm = T)
    FoldChangeLevel1vsLevel2 = ifelse(meanLevel1/meanLevel2 > 
        1, meanLevel1/meanLevel2, -meanLevel2/meanLevel1)
    output = data.frame(ID = dimnames(datExpr)[[2]], corPearson = corPearson, 
        pvalueStudent = pvalueStudent, qvalueStudent = q.Student, 
        FoldChange = FoldChangeLevel1vsLevel2, meanFirstGroup = meanLevel1, 
        meanSecondGroup = meanLevel2, AreaUnderROC = AreaUnderROC)
    if (kruskalTest) {
        output = data.frame(output, pvaluekruskal = pvaluekruskal, 
            qvaluekruskal = q.kruskal)
    }
    names(output)[-(1:4)] = paste(names(output)[-c(1:4)], levelsy[[1]], "vs", 
        levelsy[[2]], sep = ".")
    output
}

