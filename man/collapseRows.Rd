\name{collapseRows}
\alias{collapseRows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Select one representative row per group}
\description{The function collapses the rows of a numeric matrix, i.e. it selects one representative row for
each group of rows specified by the vector \code{rowGroup}. The word "collapse" reflects the fact that the
function outputs (among other things) a data frame whose rows are a subset of the original input data
\code{datET}.
The functions implements several network-based and biostatistical methods for finding a representative row
for each group specified in \code{rowGroup}. Optionally, the function identifies the representative row
according to the least number of missing data, the highest sample variance, and the highest connectivity.
One of the advantages of this function is that it implements default settings which have worked well in
numerous applications. Below, we describe these default settings in more detail. 
}
\usage{
collapseRows(datET, rowGroup, rowID, 
             method="MaxMean", connectivityBasedCollapsing=FALSE,
             methodFunction=NULL, 
             connectivityPower=1, selectFewestMissing=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datET}{matrix or data frame containing numeric values where rows correspond to variables (e.g.
microarray probes) and columns correspond to observations (e.g. microarrays). Each row of \code{datET} must
have a unique row identifier (specified in the vector \code{rowID}). The group label of each row  is encoded
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
"maxRowVariance" (default) = choose the row with the highest variance (across the columns of \code{datET}).
"MaxMean" or "MinMean" = choose the row with the highest or lowest mean value, respectively. 
"absMaxMean" or "absMinMean" = choose the row with the highest or lowest mean absolute value.
"function" = use this method for a user-input function (see the description of the argument
"methodFunction").
}
  \item{connectivityBasedCollapsing}{logical value.
If TRUE (default), groups with 3 or more corresponding rows will be represented by the row with the highest
connectivity according to a signed weighted correlation network adjacency matrix among the corresponding
rows. Recall that the connectivity is defined as the rows sum of the adjacency matrix. The signed weighted
adjacency matrix is defined as A=(0.5+0.5*COR)^power where power is determined by the argument
connectivityPower and COR denotes the matrix of pairwise Pearson correlation coefficients among the
corresponding rows.
}
  \item{methodFunction}{character string. It only needs to be specified if method="function" otherwise
its input is ignored.  Must be a function that takes a Nr x Nc matrix of numbers as input and outputs a
vector with the length Nc (e.g., colMeans).  This will then be the method used for collapsing values for
multiple rows into a single value for the row.
}
  \item{connectivityPower}{Positive number (typically integer) for specifying the threshold (power) used
to construct the signed weighted adjacency matrix, see the description of connectivityBasedCollapsing. This
option is only used if connectivityBasedCollapsing=TRUE.
}
  \item{selectFewestMissing}{ logical values. If TRUE (default), the input expression matrix is trimmed
such that for each group only the rows with the fewest number of missing values are retained.  In situations
where an equal number of values are missing (or where there is no missing data), all rows for a given group
are retained.  Whether this value is set to TRUE or FALSE, all rows with >90\% missing data are omitted from
the analysis.
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
Often several rows have the same minimum number of missing values and a representative must be chosen among
those rows.
The following remarks only apply to the case when more than 1 row have the same minimum number of missing
values.
We distinguish 2 situations:
(1) If a group corresponds to exactly 2 rows then the corresponding row with the highest variance is selected
if \code{method="maxRowVariance"}. Alternative methods can be chosen as described in the desription of the
method argument. If the two rows have the same row variance, then the row with the highest mean is selected.
In the rare case, when both rows have the same number of missing data, the same variance, and the same mean,
then the function randomly choose a representative row.
(2) If a group corresponds to more than 2 rows, then the  function calculates a signed weighted correlation
network (with power specified in connectivityPower) among the corresponding rows if
connectivityBasedCollapsing=TRUE. Next the function calculates the network connectivity of each row (closely
related to the sum or correlations with the other matching rows). Next it chooses the most highly connected
row as representative. If several probes have the same maximum connectivity, the first such row is chosen.

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
rows corresponding to the different row groups rather than individual row identifiers. }

\item{group2row}{ is a matrix whose rows correspond to the unique group labels and whose 2 columns report
which group label (first column called \code{group}) is represented by what row label (second column called
\code{selectedRowID})}.

\item{selectedRow}{ is a logical vector whose components are TRUE for probes selected as representatives and
FALSE otherwise. It has the same length as the vector probeID.}

}
\references{
Miller JA, Langfelder P, Horvath S (2010) Correlation network methods for
merging expression data from different platforms. The collapseRows function. Technical Report.
}
\author{
Jeremy A. Miller, Steve Horvath, Peter Langfelder
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

    # The code simulates a data frame (called dat1) of correlated rows.
    # You can skip this part and start at the line called Typical Input Data
    # The first column of the data frame will contain row identifiers
    # number of columns (e.g. observations or microarrays)
    m=60
    # number of rows (e.g. variables or probes on a microarray) 
    n=3000
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

} 

\keyword{misc}
	
