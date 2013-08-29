\name{branchEigengeneDissim}
\alias{branchEigengeneDissim}
\alias{mtd.branchEigengeneDissim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Branch dissimilarity based on eigennodes (eigengenes).
}
\description{
Calculation of branch dissimilarity based on eigennodes (eigengenes) in single set and multi-data
situations. This function is used as a plugin for the
dynamicTreeCut package and the user should not call this function directly. This function is experimental
and subject to change.
}
\usage{
branchEigengeneDissim(
  expr, 
  branch1, branch2, 
  corFnc = cor, corOptions = list(use = "p"), 
  signed = TRUE, ...)

mtd.branchEigengeneDissim(
  multiExpr, 
  branch1, branch2,
  corFnc = cor, corOptions = list(use = 'p'),
  consensusQuantile = 0, 
  signed = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{
Expression data.
}
  \item{multiExpr}{
Expression data in multi-set format.
}
  \item{branch1}{
Branch 1.
}
  \item{branch2}{
Branch 2.
}
  \item{corFnc}{
Correlation function.
}
  \item{corOptions}{
Other arguments to the correlation function.
}
  \item{consensusQuantile}{
Consensus quantile.
}
  \item{signed}{
Should the network be considered signed?
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
