\name{selectFewestConsensusMissing}
\alias{selectFewestConsensusMissing}
\title{
Select columns with the lowest consensus number of missing data
}
\description{
Given a \code{\link{multiData}} structure, this function calculates the consensus number of present
(non-missing) data
for each variable (column) across the data sets, forms the consensus and for each group selects variables
whose consensus proportion of present data is at least \code{selectFewestMissing} (see usage below).
}
\usage{
selectFewestConsensusMissing(
    mdx, 
    colID, 
    group, 
    minProportionPresent = 1, 
    consensusQuantile = 0, 
    verbose = 0,
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdx}{
A \code{\link{multiData}} structure. All sets must have the same columns.
}
  \item{colID}{ Character vector of column identifiers.  This must include all the column names from
\code{mdx}, but can include other values as well. Its entries must be unique (no duplicates) and no
missing values are permitted. 
}
  \item{group}{
Character vector whose components contain the group label (e.g. a character string) for
each entry of \code{colID}. This vector must be of the same length as the vector \code{colID}. In gene
expression applications, this vector could contain the gene symbol (or a co-expression module label).
}
  \item{minProportionPresent}{A numeric value between 0 and 1 (logical values will be coerced to numeric). 
Denotes the minimum consensus 
fraction of present data in each column that will result in the column being retained.
}
  \item{consensusQuantile}{A number between 0 and 1 giving the quantile probability for consensus calculation.
0 means the minimum value (true consensus) will be used.}
  \item{verbose}{
Level of verbosity; 0 means silent, larger values will cause progress messages to be printed.
}
  \item{...}{Other arguments that should be considered undocumented and subject to change.}
}
\details{
A 'consensus' of a vector (say 'x') is simply defined as the quantile with probability
\code{consensusQuantile} of the vector x. This function calculates, for each variable in \code{mdx}, its
proportion of present (i.e., non-NA and non-NaN) 
values in each of the data sets in \code{mdx}, and forms the consensus. Only
variables whose consensus proportion of present data is at least \code{selectFewestMissing} are retained.
}
\value{
A logical vector with one element per variable in \code{mdx}, giving \code{TRUE} for the retained
variables.
}
\author{
Jeremy Miller and Peter Langfelder
}
\seealso{
\code{\link{multiData}}
}
\keyword{misc }
