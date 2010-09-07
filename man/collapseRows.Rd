\name{collapseRows}
\alias{collapseRows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapse Rows in Numeric Matrix
}
\description{
This function collapses the rows of a numeric matrix from one set of variables into a less specific set.  For example, it will take as input expression data with rows as probes and output expression data with rows as genes, collapsing all probes for a given gene into one set of expression values.  Most of the parameters are used to specify which of the multiple probes for a given gene are kept.
}
\usage{
collapseRows(datET, symbolVector, probeVector, 
             method="maxRowVariance", 
             connectivityBasedCollapsing=TRUE,
             methodFunction=NULL, 
             connectivityPower=1, 
             selectFewestMissing=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datET}{
Matrix (or data frame) of numeric values (typically expression data) where rows are probes and columns are samples.
}
  \item{symbolVector}{
A character vector of gene symbols which correspond to the probes in "probeVector"
}
  \item{probeVector}{
A character vector of probes.  This should include all the probes from rownames(datET), but can include other probes.
}
  \item{method}{
Built-in options for possible methods for determining the "best" probe to keep.  These are the options:
"maxRowVariance" (default) = choose the probe with the highest variance;
"Max" / "Min" = choose the probe with the highest / lowest mean expression level
"absMax" / "absMin" = choose the probe with the highest / lowest mean absolute value of the expression level
function = use this method for a user-input function (see \code{methodFunction} parameter)
}
  \item{connectivityBasedCollapsing}{
If TRUE (default), genes with 3+ corresponding probes will be collapsed such that values corresponding to the probe with the highest connectivity (signed adjacency raised to a power of connectivityPower) will be used.
}
  \item{methodFunction}{
Should only be set if \code{method="function"}.  Must be a function that takes a Nr x Nc matrix of numbers
as input and outputs a vector with the length Nc (ie, colMeans).  This will then be the method used for
collapsing values for multiple probes into a single value for the gene. 
}
  \item{connectivityPower}{
Only used if connectivityBasedCollapsing=TRUE (see connectivityBasedCollapsing parameter).  Value must be a number > 0.
}
  \item{selectFewestMissing}{
If TRUE (default), the input expression matrix is trimmed such that only the probes for each gene with the
fewest number of missing probes are retained.  In situations where an equal number of probes are missing (or
where there is no missing data), all probes for a given gene are retained.  Whether this value is set to
TRUE or FALSE, all probes with >90\% missing data are omitted from the analysis.
}

}
\details{
The function is robust to missing data.  Also, if probe names are missing, they are inferred where possible.
}
\value{
The output is an expression matrix with the same columns as datET, but with rows corresponding to genes rather than probes.
%% ...
}
\author{
Jeremy A. Miller
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\keyword{ misc }
