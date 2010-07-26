\name{collapseRows}
\alias{collapseRows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapse Rows in Numeric Matrix
}
\description{
This function collapses the rows of a numeric matrix from one set of variables into a less specific set.  For example, it will take as input expression data with rows as probes and output expression data with rows as genes, collapsing all probes for a given gene into one set of expression values.
}
\usage{
collapseRows(datET, symbolVector, probeVector, method = max, connectivityPower=NULL)
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
Must be a function that takes a vector of numbers as input and outputs a single number (ie, min, max, mean, median).  For genes with 2 probes (also 3+ if connectivityPower=NULL), this is the method used for collapsing values for multiple probes into a single value for the gene.
}
  \item{connectivityPower}{
If a value is inputted, genes with 3+ corresponding probes will be collapsed such that values corresponding to the probe with the highest connectivity (unsigned adjacency raised to a power of connectivityPower) will be used.  Value must be a number > 0.
}

}
\details{
The function is robust to missing data, in that NAs are only outputted if expression values for ALL probes for a given gene are missing.  Also, if probe names are missing, they are inferred where possible.
}
\value{
The output is an expression matrix with the same columns as datET, but with rows corresponding to genes rather than probes.
}
\author{
Jeremy A. Miller
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\keyword{ misc }
