\name{branchSplit}
\alias{branchSplit}

\title{
Branch split.
}
\description{
Calculation of branch split based on expression data. This function is used as a plugin for the
dynamicTreeCut package and the user should not call this function directly.
}
\usage{
branchSplit(
  expr, 
  branch1, branch2, 
  discardProp = 0.05, minCentralProp = 0.75, 
  nConsideredPCs = 3, 
  signed = FALSE, 
  getDetails = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{
Expression data.
}
  \item{branch1}{
Branch 1,
}
  \item{branch2}{
Branch 2.
}
  \item{discardProp}{
Proportion of data to be discarded as outliers.
}
  \item{minCentralProp}{
Minimum central proportion
}
  \item{nConsideredPCs}{
Number of principal components to consider.
}
  \item{signed}{
Should the network be considered signed?
}
  \item{getDetails}{
Should details of the calculation be returned?
}
  \item{\dots}{
Other arguments. Present for compatibility; currently unusued.
}
}
\value{
A single number or a list containing detils of the calculation.
}
\author{
Peter Langfelder
}
\keyword{ misc }% __ONLY ONE__ keyword per line
