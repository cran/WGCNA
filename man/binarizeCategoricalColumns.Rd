\name{binarizeCategoricalColumns}
\alias{binarizeCategoricalColumns}
\alias{binarizeCategoricalColumns.pairwise}
\alias{binarizeCategoricalColumns.forRegression}
\alias{binarizeCategoricalColumns.forPlots}
\title{
Turn categorical columns into sets of binary indicators 
}
\description{
Given a data frame with (some) categorical columns, 
this function creates a set of indicator variables for the various possible
sets of levels.
}
\usage{
binarizeCategoricalColumns(
   data,
   convertColumns = NULL,
   considerColumns = NULL,
   maxOrdinalLevels = 3,
   levelOrder = NULL,
   minCount = 3,
   val1 = 0, val2 = 1,
   includePairwise = FALSE,
   includeLevelVsAll = TRUE,
   dropFirstLevelVsAll = TRUE,
   dropUninformative = TRUE,
   includePrefix = TRUE,
   prefixSep = ".",
   nameForAll = "all",
   levelSep = NULL,
   levelSep.pairwise = if (length(levelSep)==0) ".vs." else levelSep,
   levelSep.vsAll = if (length(levelSep)==0) 
              (if (nameForAll=="") "" else ".vs.") else levelSep,
   checkNames = FALSE,
   includeLevelInformation = FALSE)

binarizeCategoricalColumns.pairwise(
   data, 
   maxOrdinalLevels = 3,
   convertColumns = NULL,
   considerColumns = NULL,
   levelOrder = NULL,
   val1 = 0, val2 = 1, 
   includePrefix = TRUE,
   prefixSep = ".", 
   levelSep = ".vs.",
   checkNames = FALSE)

binarizeCategoricalColumns.forRegression(
   data, 
   maxOrdinalLevels = 3,
   convertColumns = NULL,
   considerColumns = NULL,
   levelOrder = NULL,
   val1 = 0, val2 = 1,
   includePrefix = TRUE,
   prefixSep = ".",
   checkNames = TRUE)

binarizeCategoricalColumns.forPlots(
   data, 
   maxOrdinalLevels = 3,
   convertColumns = NULL,
   considerColumns = NULL,
   levelOrder = NULL,
   val1 = 0, val2 = 1,
   includePrefix = TRUE,
   prefixSep = ".",
   checkNames = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame.
}
  \item{convertColumns}{
Optional character vector giving the column names of the columns to be converted. See \code{maxOrdinalLevels}
below.
}
  \item{considerColumns}{
Optional character vector giving the column names of columns that should be looked at and possibly converted.
If not given, all columns will be considered. See \code{maxOrdinalLevels} below. 
}
  \item{maxOrdinalLevels}{
When \code{convertColumns} above is \code{NULL}, the function looks at all columns in \code{considerColumns}
and converts all non-numeric columns and those numeric columns that have at most \code{maxOrdinalLevels}
unique values. A column is considered numeric if its storage mode is numeric or if it is character and all
entries with the expception of "NA", "NULL" and "NO DATA" represent valid numbers.
}
  \item{levelOrder}{
Optional list giving the ordering of levels (unique values) in each of the converted columns. Best used in
conjunction with \code{convertColumns}.
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
  \item{includePrefix}{
Logical: should the column name of the binarized column be included in column names of the output? See
details.
}
  \item{prefixSep}{
Separator of column names and level names in column names of the output. See details.
}
  \item{nameForAll}{ Character string that represents "all others" in the column names of indicators of level vs. all
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

\code{binarizeCategoricalColumns} is the most general function, the rest are convenience wrappers that set some of the
options to achieve the following:

\code{binarizeCategoricalColumns.pairwise} returns only pairwise (level vs. level) binary indicators.

\code{binarizeCategoricalColumns.forRegression} returns only level vs. all others binary indicators, with the first
(according to \code{levelOrder}) level
vs. all removed. This is essentially the same as would be returned by \code{\link{model.matrix}} except for the column
representing intercept.

\code{binarizeCategoricalColumns.forPlots} returns only level vs. all others binary indicators and keeps them all.

The columns to be converted are identified as follows. If \code{considerColumns} is given, columns not
contained in it will not be converted, even if they are included in \code{convertColumns}. 

If \code{convertColumns} is given, those columns will
be converted (except any not contained in non-empty \code{considerColumns}). If \code{convertColumns} 
is \code{NULL}, the function converts columns that are not numeric (as reported by \code{\link{is.numeric}}) and those
numeric columns that have at most \code{maxOrdinalValues} unique non-missing values. 

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

Internally, the function calls \code{\link{binarizeCategoricalVariable}} for each column that is converted.

}
\value{
A data frame in which the converted columns have been replaced by sets of binarized indicators. When
\code{includeLevelInformation} is
\code{TRUE}, the attribute \code{includedLevels} is a table with one column per output column and two rows,
giving the two levels (unique values of x) represented by the column.


}
\author{
Peter Langfelder
}
\examples{
set.seed(2);
x = data.frame(a = sample(c("A", "B", "C"), 15, replace = TRUE),
               b = sample(c(1:3), 15, replace = TRUE));
out = binarizeCategoricalColumns(x, includePairwise = TRUE, includeLevelVsAll = TRUE,
                     includeLevelInformation = TRUE);
data.frame(x, out);
attr(out, "includedLevels")

}
\keyword{misc}
