\name{binarizeCategoricalVariable}
\alias{binarizeCategoricalVariable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turn a categorical variable into a set of binary indicators
}
\description{
Given a categorical variable, this function creates a set of indicator variables for the various possible
sets of levels.
}
\usage{
binarizeCategoricalVariable(
   x,
   levelOrder = NULL,
   ignore = NULL,
   minCount = 3,
   val1 = 0, val2 = 1,
   includePairwise = TRUE,
   includeLevelVsAll = FALSE,
   dropFirstLevelVsAll = FALSE,
   dropUninformative = TRUE,
   namePrefix = "",
   levelSep = NULL,
   nameForAll = "all",
   levelSep.pairwise = if (length(levelSep)==0) ".vs." else levelSep,
   levelSep.vsAll = if (length(levelSep)==0) 
                       (if (nameForAll=="") "" else ".vs.") else levelSep,
   checkNames = FALSE,
   includeLevelInformation = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector with categorical values.
}
  \item{levelOrder}{
Optional specification of the levels (unique values) of \code{x}. Defaults to sorted unique values of
\code{x}, but can be used to only include a subset of the existing levels as well as to specify the order of
the levels in the output variables.
}
  \item{ignore}{
Optional specification of levels of \code{x} that are to be ignored. Note that the levels are ignored only when
deciding which variables to include in the output; the samples with these values of \code{x} will be included
in "all" in indicators of level vs. all others.
}
  \item{minCount}{
Levels of \code{x} for which there are fewer than \code{minCount} elements will be ignored.
}
  \item{val1}{
Value for the lower level in binary comparisons.
}
  \item{val2}{
Value for the higher level in binary comparisons.
}
  \item{includePairwise}{
Logical: should pairwise binary indicators be included? For each pair of levels, the indicator is \code{val1}
for the lower level (earlier in \code{levelOrder}), \code{val2} for the higher level and \code{NA} otherwise.
}
  \item{includeLevelVsAll}{
Logical: should binary indicators for each level be included? The indicator is \code{val2} where \code{x}
equals the level and \code{val1} otherwise. 
}
  \item{dropFirstLevelVsAll}{
Logical: should the column representing first level vs. all be dropped? This makes the resulting matrix of
indicators usable for regression models. 
}
  \item{dropUninformative}{
Logical: should uninformative (constant) columns be dropped?
}
  \item{namePrefix}{
Prefix to be used in column names of the output.
}
\item{nameForAll}{
When naming columns that represent a level vs. all others, \code{nameForAll} will be used to represent all
others.
}
  \item{levelSep}{
Separator for levels to be used in column names of the output. If \code{NULL}, pairwise and level vs. all indicators will
use different level separators set by \code{levelSep.pairwise} and \code{levelSep.vsAll}.
}
  \item{levelSep.pairwise}{
Separator for levels to be used in column names for pairwise indicators in the output.
}
  \item{levelSep.vsAll}{
Separator for levels to be used in column names for level vs. all indicators in the output.
}
  \item{checkNames}{
Logical: should the names of the output be made into syntactically correct R language names?
}
  \item{includeLevelInformation}{
Logical: should information about which levels are represented by which columns be included in the attributes
of the output?
}
}
\details{
The function creates two types of indicators. The first is one level (unique value) of \code{x} vs. all
others, i.e., for a given level, the indicator is \code{val2} (usually 1) for all elements of \code{x} that
equal the level, and \code{val1} (usually 0) 
otherwise. Column names for these indicators are the concatenation of \code{namePrefix}, the level,
\code{nameSep} and \code{nameForAll}. The level vs. all indicators are created for all levels that have at
least \code{minCounts} samples, are present in \code{levelOrder} (if it is non-NULL) and are not included in
\code{ignore}.

The second type of indicator encodes binary comparisons. For each pair of levels (both with at least
\code{minCount} samples), the indicator is \code{val2} (usually 1) for the higher level and \code{val1}
(usually 0) for the lower level. The level order is given by \code{levelOrder} (which defaults to the sorted
levels of \code{x}), assumed to be sorted in increasing order. All levels with at least \code{minCount}
samples that are included in \code{levelOrder} and not included in \code{ignore} are included.

}
\value{
A matrix containing the indicators variabels, one in each column. When \code{includeLevelInformation} is
\code{TRUE}, the attribute \code{includedLevels} is a table with one column per output column and two rows,
giving the two levels (unique values of x) represented by the column.
}
\author{
Peter Langfelder
}
\seealso{
Variations and wrappers for this function:
\code{binarizeCategoricalColumns} for binarizing several columns of a matrix or data frame
}
\examples{
set.seed(2);
x = sample(c("A", "B", "C"), 15, replace = TRUE);
out = binarizeCategoricalVariable(x, includePairwise = TRUE, includeLevelVsAll = TRUE);
data.frame(x, out);
attr(out, "includedLevels")
# A different naming for level vs. all columns
binarizeCategoricalVariable(x, includeLevelVsAll = TRUE, nameForAll = "");
}
\keyword{misc}
