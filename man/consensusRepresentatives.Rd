\name{consensusRepresentatives}
\alias{consensusRepresentatives}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Consensus selection of group representatives
}
\description{
Given multiple data sets corresponding to the same variables and a grouping of variables into groups, 
the function selects a
representative variable for each group using a variety of possible selection approaches. Typical uses include
selecting a representative probe for each gene in microarray data.
}
\usage{
consensusRepresentatives(
   mdx, 
   group, 
   colID, 
   consensusQuantile = 0, 
   method = "MaxMean", 
   useGroupHubs = TRUE, 
   calibration = c("none", "full quantile"), 
   selectionStatisticFnc = NULL, 
   connectivityPower = 1, 
   minProportionPresent = 1, 
   getRepresentativeData = TRUE, 
   statisticFncArguments = list(), 
   adjacencyArguments = list(), 
   verbose = 2, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdx}{ A \code{\link{multiData}} structure. All sets must have the same columns.
}

  \item{group}{
Character vector whose components contain the group label (e.g. a character string) for
each entry of \code{colID}. This vector must be of the same length as the vector \code{colID}. In gene
expression applications, this vector could contain the gene symbol (or a co-expression module label).
}

  \item{colID}{Character vector of column identifiers.  This must include all the column names from
\code{mdx}, but can include other values as well. Its entries must be unique (no duplicates) and no
missing values are permitted.
}

  \item{consensusQuantile}{
A number between 0 and 1 giving the quantile probability for consensus calculation.
0 means the minimum value (true consensus) will be used.}

  \item{method}{character string for determining which method is used to choose the representative 
(when \code{useGroupHubs} is \code{TRUE}, this method is only used for groups with 2
variables).
The following values can be used:
"MaxMean" (default) or "MinMean" return the variable with the highest or lowest mean value, respectively;
"maxRowVariance" return the variable with the highest variance;
"absMaxMean" or "absMinMean" return the variable with the highest or lowest mean absolute value; and
"function" will call a user-input function (see the description of the argument
\code{selectionStatisticFnc}). The built-in functions can be instructed to use robust analogs (median and
median absolute deviation) by also specifying \code{statisticFncArguments=list(robust = TRUE)}.
}
  \item{useGroupHubs}{Logical: if \code{TRUE}, groups with 3 or more variables will be
represented by the variable with the highest
connectivity according to a signed weighted correlation network adjacency matrix among the corresponding
rows. The connectivity is defined as the row sum of the adjacency matrix. The signed weighted
adjacency matrix is defined as A=(0.5+0.5*COR)^power where power is determined by the argument
\code{connectivityPower} and COR denotes the matrix of pairwise correlation coefficients among the
corresponding rows. Additional arguments to the underlying function \code{\link{adjacency}} can be specified
using the argument \code{adjacencyArguments} below.
}
  \item{calibration}{Character string describing the method of calibration of the selection statistic among
the data sets. Recognized values are \code{"none"} (no calibration) and \code{"full quantile"} (quantile
normalization). }
  \item{selectionStatisticFnc}{User-supplied function used to calculate the selection statistic when
\code{method} above equals \code{"function"}.  The function must take argumens \code{x} (a matrix) and
possibly other arguments that can be specified using \code{statisticFncArguments} below. The return value
must be a vector with one component per column of \code{x} giving the selection statistic for each column. 
}
  \item{connectivityPower}{Positive number (typically integer) for specifying the soft-thresholding power used
to construct the signed weighted adjacency matrix, see the description of \code{useGroupHubs}.
This option is only used if \code{useGroupHubs} is \code{TRUE}.
}
  \item{minProportionPresent}{
A number between 0 and 1 specifying a filter of candidate probes. Specifically, for each group, the variable 
with the maximum consensus proportion of present data is found. Only variables whose consensus proportion of 
present data is at least \code{minProportionPresent} times the  maximum consensus proportion are retained as 
candidates for being a representative.
}
  \item{getRepresentativeData}{Logical: should the representative data, i.e., \code{mdx} restricted to
the representative variables, be returned? }
  \item{statisticFncArguments}{ A list giving further arguments to the selection statistic function. Can be
used to supply additional arguments to the user-specified \code{selectionStatisticFnc}; the value 
\code{list(robust = TRUE)} can be used with the built-in functions to use their robust variants.}
  \item{adjacencyArguments}{Further arguments to the function \code{adjacency}, e.g.
\code{adjacencyArguments=list(corFnc = "bicor", corOptions = "use = 'p', maxPOutliers = 0.05")} will select
the robust correlation \code{\link{bicor}} with a good set of options. Note that the \code{\link{adjacency}}
arguments \code{type} and \code{power} cannot be changed.
}
  \item{verbose}{
Level of verbosity; 0 means silent, larger values will cause progress messages to be printed.
}
  \item{indent}{
Indent for the diagnostic messages; each unit equals two spaces.
}
}
\details{
This function was inspired by \code{\link{collapseRows}}, but there are also important differences. This function
focuses on selecting representatives; when summarization is more important, \code{collapseRows} provides more
flexibility since it does not require that a single representative be selected. 

This function and \code{collapseRows} use different input and ouput conventions; user-specified functions need
to be tailored differently for \code{collapseRows} than for \code{consensusRepresentatives}.

Missing data are allowed and are treated as missing at random. If \code{rowID} is \code{NULL}, it is replaced
by the variable names in \code{mdx}.

All groups with a single variable are represented by that variable, unless the consensus proportion of present
data in the variable is lower than \code{minProportionPresent}, in which case the variable and the group are
excluded from the output. 

For all variables belonging to groups with 2 variables (when \code{useGroupHubs=TRUE}) or with at least 2 variables
(when \code{useGroupHubs=FALSE}), selection statistics are calculated in each set (e.g., the selection
statistic may be the mean, variance, etc). This results in a matrix of selection statistics (one entry per
variable per data set). The selection statistics are next optionally calibrated (normalized) between sets to
make them comparable; currently the only implemented calibration method is quantile normalization. 

For
each variable, the consensus selection statistic is defined as the 
consensus of the (calibrated) selection statistics across the data sets is calculated. The
'consensus' of a vector (say 'x') is simply defined as the quantile with probability
\code{consensusQuantile} of the vector x. Important exception: for the \code{"MinMean"} and
\code{"absMinMean"} methods, the consensus is the quantile with probability \code{1-consensusQuantile}, since
the idea of the consensus is to select the worst (or close to worst) value across the data sets.

For each group, the representative is selected as the variable with the best (typically highest, but for
\code{"MinMean"} and 
\code{"absMinMean"} methods the lowest) consensus selection statistic.

If \code{useGroupHubs=TRUE}, the intra-group connectivity is calculated for all variables in each set. The
intra-group connectivities are optionally calibrated (normalized) between sets, and consensus intra-group
connectivity is calculated similarly to the consensus selection statistic above. In each group, the variable
with the highest consensus intra-group connectivity is chosen as the representative. 

}
\value{
\item{representatives}{A named vector giving, for each group, the selected representative (input \code{rowID}
or the variable (column) name in \code{mdx}). Names correspond to groups.}
\item{varSelected}{A logical vector with one entry per variable (column) in input \code{mdx} (possibly
after restriction to variables occurring in \code{colID}), \code{TRUE} if the column was selected as a
representative.}
\item{representativeData}{Only present if \code{getRepresentativeData} is \code{TRUE};
the input \code{mdx} restricted to the representative variables, with column
names changed to the corresponding groups.}
}
\author{
Peter Langfelder, based on code by Jeremy Miller
}
\seealso{
\code{\link{multiData}} for a description of the \code{multiData} structures;
\code{\link{collapseRows}} that solves a related but different problem. Please note the differences in input
and output!
}
\keyword{misc}
