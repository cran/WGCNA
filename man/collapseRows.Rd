\name{collapseRows}
\alias{collapseRows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Select one representative row per group}
\description{Abstractly speaking, the function allows one to collapse the rows of a numeric matrix, 
e.g. by forming an average or selecting one representative row for each group of rows specified by a 
grouping variable (referred to as \code{rowGroup}). The word "collapse" reflects the fact that the 
method yields a new matrix whose rows correspond to other rows of the original input data. The function 
implements several network-based and biostatistical methods for finding a representative row for each 
group specified in \code{rowGroup}.
Optionally, the function identifies the representative row according to the least number of missing data, 
the highest sample mean, the highest sample variance, the highest connectivity. One of the advantages of 
this function is that it implements default settings which have worked well in numerous applications. 
Below, we describe these default settings in more detail. 
}
\usage{
collapseRows(datET, rowGroup, rowID, 
             method="MaxMean", connectivityBasedCollapsing=FALSE,
             methodFunction=NULL, connectivityPower=1,
             selectFewestMissing=TRUE, thresholdCombine=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datET}{matrix or data frame containing numeric values where rows correspond to variables (e.g.
microarray probes) and columns correspond to observations (e.g. microarrays). Each row of \code{datET} must
have a unique row identifier (specified in the vector \code{rowID}). The group label of each row is encoded
in the vector \code{rowGroup}. While \code{rowID} should have non-missing, unique values (identifiers), the
values of the vector \code{rowGroup} will typically not be unique since the function aims to pick a
representative row for each group.
}
  \item{rowGroup}{character vector whose components contain the group label (e.g. a character string) for
each row of \code{datET}. This vector needs to have the same length as the vector \code{rowID}. In gene
expression applications, this vector could contain the gene symbol (or a co-expression module label).
}
  \item{rowID}{character vector of row identifiers.  This should include all the rows from
rownames(\code{datET}), but can include other rows. Its entries should be unique (no duplicates) and no
missing values are permitted. If the row identifier is missing for a given row, we suggest you remove this
row from \code{datET} before applying the function.
}
\item{method}{character string for determining which method is used to choose a probe among
exactly 2 corresponding rows or when connectivityBasedCollapsing=FALSE. These are the options:
"MaxMean" (default) or "MinMean" = choose the row with the highest or lowest mean value, respectively.
"maxRowVariance" = choose the row with the highest variance (across the columns of \code{datET}). 
"absMaxMean" or "absMinMean" = choose the row with the highest or lowest mean absolute value.
"ME" = choose the eigenrow (first principal component of the rows in each group).  Note that with this
       method option, connectivityBasedCollapsing is automatically set to FALSE.
"Average" = for each column, take the average value of the rows in each group 
"function" = use this method for a user-input function (see the description of the argument
"methodFunction").
Note: if method="ME", "Average" or "function", the output parameters "group2row" and "selectedRow" are not informative.
}
  \item{connectivityBasedCollapsing}{logical value.
If TRUE, groups with 3 or more corresponding rows will be represented by the row with the highest
connectivity according to a signed weighted correlation network adjacency matrix among the corresponding
rows. Recall that the connectivity is defined as the rows sum of the adjacency matrix. The signed weighted
adjacency matrix is defined as A=(0.5+0.5*COR)^power where power is determined by the argument
\code{connectivityPower} and COR denotes the matrix of pairwise Pearson correlation coefficients among the
corresponding rows.
}
  \item{methodFunction}{character string. It only needs to be specified if method="function" otherwise
its input is ignored.  Must be a function that takes a Nr x Nc matrix of numbers as input and outputs a
vector with the length Nc (e.g., colMeans).  This will then be the method used for collapsing values for
multiple rows into a single value for the row.
}
  \item{connectivityPower}{Positive number (typically integer) for specifying the threshold (power) used
to construct the signed weighted adjacency matrix, see the description of \code{connectivityBasedCollapsing}.
This option is only used if connectivityBasedCollapsing=TRUE.
}
  \item{selectFewestMissing}{logical values. If TRUE (default), the input expression matrix is trimmed
such that for each group only the rows with the fewest number of missing values are retained.  In situations
where an equal number of values are missing (or where there is no missing data), all rows for a given group
are retained.  Whether this value is set to TRUE or FALSE, all rows with >90\% missing data are omitted from
the analysis.
}
  \item{thresholdCombine}{Number between -1 and 1, or NA. If NA (default), this input is ignored.
If a number between -1 and 1 is input, this value is taken as a threshold value, and collapseRows proceeds
following the "maxMean" method, but ONLY for ids with correlations of R>thresholdCombine.  Specifically:
...1) If there is one id/group, keep the id
...2) If there are 2 ids/group, take the maximum mean expression if their correlation is > thresholdCombine
...3) If there are 3+ ids/group, iteratively repeat (2) for the 2 ids with the highest correlation until 
      all ids remaining have correlation < thresholdCombine for each group
Note that this option usually results in more than one id per group; therefore, one must use care when 
implementing this option for use in comparisons between multiple matrices / data frames.
}

}
\details{
The function is robust to missing data.  Also, if rowIDs are missing, they are inferred according to the
rownames of datET when possible.  
When a group corresponds to only 1 row then it is represented by this row since there is no other choice.
Having said this, the row may be removed if it contains an excessive amount of missing data (90 percent or
more missing values), see the description of the argument \code{selectFewestMissing} for more details.

A group is represented by a corresponding row with the fewest number of missing data if
\code{selectFewestMissing} has been set to TRUE. 
Often several rows have the same minimum number of missing values (or no missing values) and a representative 
must be chosen among those rows.  In this case we distinguish 2 situations:
(1) If a group corresponds to exactly 2 rows then the corresponding row with the highest average is selected
if \code{method="maxMean"}. Alternative methods can be chosen as described in \code{method}.
(2) If a group corresponds to more than 2 rows, then the function calculates a signed weighted correlation
network (with power specified in \code{connectivityPower}) among the corresponding rows if
\code{connectivityBasedCollapsing=TRUE}. Next the function calculates the network connectivity of each row 
(closely related to the sum or correlations with the other matching rows). Next it chooses the most highly
connected row as representative. If connectivityBasedCollapsing=FALSE, then \code{method} is used. 
For both situations, if more than one row has the same value, the first such row is chosen.

Setting \code{thresholdCombine} is a special case of this function, as not all ids for a single group are
necessarily collapsed--only those with similar expression patterns are collapsed.  We suggest
using this option when the goal is to decrease the number of ids for computational reasons, but when
ALL ids for a single group should not be combined (for example, if two probes could represent different 
splice variants for the same gene for many genes on a microarray).

Example application: when dealing with microarray gene expression data then the rows of \code{datET} may
correspond to unique probe identifiers and \code{rowGroup} may contain corresponding gene symbols. Recall
that multiple probes (specified using \code{rowID}=ProbeID) may correspond to the same gene symbol
(specified using \code{rowGroup}=GeneSymbol). In this case, \code{datET} contains the input expression data
with rows as rowIDs and output expression data with rows as gene symbols, collapsing all probes for a given
gene symbol into one representative.  

}
\value{
The output is a list with the following components.

\item{datETcollapsed}{ is a numeric matrix with the same columns as the input matrix \code{datET}, but with
rows corresponding to the different row groups rather than individual row identifiers.  (If thresholdCombine
is set, then rows still correspond to individual row identifiers.) }

\item{group2row}{ is a matrix whose rows correspond to the unique group labels and whose 2 columns report
which group label (first column called \code{group}) is represented by what row label (second column called
\code{selectedRowID}).  Set to NULL if method="ME" or "function".}.

\item{selectedRow}{ is a logical vector whose components are TRUE for probes selected as representatives and
FALSE otherwise. It has the same length as the vector probeID.  Set to NULL if method="ME" or "function".}

}
\references{
Miller JA, Langfelder P, Cai C, Horvath S (2010) Strategies for optimally aggregating gene expression data: 
The collapseRows R function. Technical Report.
}
\author{
Jeremy A. Miller, Steve Horvath, Peter Langfelder, Chaochao Cai
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
    ########################################################################
    # EXAMPLE 1:
    # The code simulates a data frame (called dat1) of correlated rows.
    # You can skip this part and start at the line called Typical Input Data
    # The first column of the data frame will contain row identifiers
    # number of columns (e.g. observations or microarrays)
    m=60
    # number of rows (e.g. variables or probes on a microarray) 
    n=500
    # seed module eigenvector for the simulateModule function
    MEtrue=rnorm(m)
    # numeric data frame of n rows and m columns
    datNumeric=data.frame(t(simulateModule(MEtrue,n)))
    RowIdentifier=paste("Probe", 1:n, sep="")
    ColumnName=paste("Sample",1:m, sep="")
    dimnames(datNumeric)[[2]]=ColumnName
    # Let us now generate a data frame whose first column contains the rowID
    dat1=data.frame(RowIdentifier, datNumeric)
    #we simulate a vector with n/5 group labels, i.e. each row group corresponds to 5 rows
    rowGroup=rep(  paste("Group",1:(n/5),  sep=""), 5 )
    
    # Typical Input Data 
    # Since the first column of dat1 contains the RowIdentifier, we use the following code
    datET=dat1[,-1]
    rowID=dat1[,1]
    
    # assign row names according to the RowIdentifier 
    dimnames(datET)[[1]]=rowID
    # run the function and save it in an object
    
    collapse.object=collapseRows(datET=datET, rowGroup=rowGroup, rowID=rowID)
    
    # this creates the collapsed data where 
    # the first column contains the group name
    # the second column reports the corresponding selected row name (the representative)
    # and the remaining columns report the values of the representative row
    dat1Collapsed=data.frame( collapse.object$group2row, collapse.object$datETcollapsed)
    dat1Collapsed[1:5,1:5]

    ########################################################################
    # EXAMPLE 2:
    # Using the same data frame as above, run collapseRows with a user-inputted function.
    # In this case we will use the mean.  Note that since we are choosing some combination
    #   of the probe values for each gene, the group2row and selectedRow output 
    #   parameters are not meaningful.

    collapse.object.mean=collapseRows(datET=datET, rowGroup=rowGroup, rowID=rowID, 
          method="function", methodFunction=colMeans)[[1]]

    # Note that in this situation, running the following code produces the identical results:

    collapse.object.mean.2=collapseRows(datET=datET, rowGroup=rowGroup, rowID=rowID,
          method="Average")[[1]]

    ########################################################################
    # EXAMPLE 3:
    # Using collapseRows to calculate the module eigengene.
    # First we create some sample data as in example 1 (or use your own!)
    m=60
    n=500
    MEtrue=rnorm(m)
    datNumeric=data.frame(t(simulateModule(MEtrue,n)))

    # In this example, rows are genes, and groups are modules.
    RowIdentifier=paste("Gene", 1:n, sep="")
    ColumnName=paste("Sample",1:m, sep="")
    dimnames(datNumeric)[[2]]=ColumnName
    dat1=data.frame(RowIdentifier, datNumeric)
    # We simulate a vector with n/100 modules, i.e. each row group corresponds to 100 rows
    rowGroup=rep(  paste("Module",1:(n/100),  sep=""), 100 )
    datET=dat1[,-1]
    rowID=dat1[,1]
    dimnames(datET)[[1]]=rowID

    # run the function and save it in an object
    collapse.object.ME=collapseRows(datET=datET, rowGroup=rowGroup, rowID=rowID, method="ME")[[1]]
    
    # Note that in this situation, running the following code produces the identical results:
    collapse.object.ME.2 = t(moduleEigengenes(expr=t(datET),colors=rowGroup)$eigengene)
    colnames(collapse.object.ME.2) = ColumnName
    rownames(collapse.object.ME.2) = sort(unique(rowGroup))
} 

\keyword{misc}
