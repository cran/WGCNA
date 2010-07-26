\name{keepCommonProbes}
\alias{keepCommonProbes}
\title{ Keep probes that are shared among given data sets }
\description{
  This function strips out probes that are not shared by all given data sets, and orders the remaining
common probes using the same order in all sets.
}
\usage{
keepCommonProbes(multiExpr, orderBy = 1)
}
\arguments{
  \item{multiExpr}{  expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{orderBy}{ index of the set by which probes are to be ordered. }
}
\value{
  Expression data in the same format as the input data, containing only common probes.
  
}
\author{ Peter Langfelder }
\seealso{\code{\link{checkSets}}}

\keyword{ misc }% __ONLY ONE__ keyword per line
