\name{mtd.subset}
\alias{mtd.subset}
\title{
Subset rows and columns in a multiData structure
}
\description{
The function restricts each \code{data} component to the given columns and rows.
}
\usage{
mtd.subset(
  # Input
  multiData, 

  # Rows and columns to keep
  rowIndex = NULL, colIndex = NULL, 
  invert = FALSE,

  # Strict or permissive checking of structure?
  permissive = FALSE, 

  # Output formatting options
  drop = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiData}{
A multiData structure.
}
  \item{rowIndex}{
A list in which each component corresponds to a set and is a vector giving the rows to be retained in that
set. All indexing methods recognized by R can be used (numeric,
logical, negative indexing, etc). If \code{NULL}, all columns will be retained in each set. Note that
setting individual elements of \code{rowIndex} to \code{NULL} will lead to errors. 
}
  \item{colIndex}{
A vector giving the columns to be retained. All indexing methods recognized by R can be used (numeric,
logical, negative indexing, etc). In addition, column names of the retained columns may be given; if a given
name cannot be matched to a column, an error will be thrown. If \code{NULL}, all columns will be retained.
}

  \item{invert}{Logical: should the selection be inverted?}

  \item{permissive}{Logical: should the function tolerate "loose" \code{multiData} input? Note that the subsetting
may lead to cryptic errors if the input \code{multiData} does not follow the "strict" format. }

  \item{drop}{Logical: should dimensions with extent 1 be dropped? }
}
\details{
A multiData structure is intended to store (the same type of) data for multiple, possibly independent,
realizations
(for example, expression data for several independent experiments). It is a list where
each component corresponds to an (independent) data set. Each component is in turn a list that can hold
various types of information but must have a \code{data} component. In a "strict" multiData structure, the
\code{data} components are required to each be a matrix or a data frame and have the same number of
columns. In a "loose" multiData structure, the \code{data} components can be anything (but for most
purposes should be of comparable type and content).

This function assumes a "strict" multiData structure unless \code{permissive} is \code{TRUE}.
}
\value{
A multiData structure containing the selected rows and columns. Attributes (except possibly dimensions and the
corresponding dimnames) are retained.
}
\author{
Peter Langfelder
}

\seealso{
 \code{\link{multiData}} to create a multiData structure.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{ misc}
