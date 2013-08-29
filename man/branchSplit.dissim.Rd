\name{branchSplit.dissim}
\alias{branchSplit.dissim}
\title{
Branch split based on dissimilarity.
}
\description{
Calculation of branch split based on a dissimilarity matrix. This function is used as a plugin for the
dynamicTreeCut package and the user should not call this function directly. This function is experimental
and subject to change.
}
\usage{
branchSplit.dissim(
  dissimMat, 
  branch1, branch2, 
  upperP, 
  minNumberInSplit = 5, 
  getDetails = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dissimMat}{
Dissimilarity matrix.}
  \item{branch1}{
Branch 1.
}
  \item{branch2}{
Branch 2.
}
  \item{upperP}{
Percentile of (closest) objects to be considered.
}
  \item{minNumberInSplit}{
Minimum number of objects to be considered.
}
  \item{getDetails}{
Should details of the calculation be returned?
}
  \item{\dots}{
Other arguments for compatibility; currently unused.
}
}
\value{
A single number or a list containing details of the calculation.
}
\author{
Peter Langfelder
}
\keyword{ misc }% __ONLY ONE__ keyword per line
