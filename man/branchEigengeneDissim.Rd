\name{branchEigengeneDissim}
\alias{branchEigengeneDissim}
\alias{branchEigengeneSimilarity}
\alias{mtd.branchEigengeneDissim}
\alias{hierarchicalBranchEigengeneDissim}
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

branchEigengeneSimilarity(
  expr, 
  branch1, 
  branch2, 
  networkOptions, 
  returnDissim = TRUE, ...)

mtd.branchEigengeneDissim(
  multiExpr, 
  branch1, branch2,
  corFnc = cor, corOptions = list(use = 'p'),
  consensusQuantile = 0, 
  signed = TRUE, reproduceQuantileError = FALSE, ...)

hierarchicalBranchEigengeneDissim(
    multiExpr,
    branch1, branch2,
    networkOptions,
    consensusTree, ...)
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

  \item{reproduceQuantileError}{Logical: should an error in the calculation from previous versions, which
caused the true consensus quantile to be \code{1-consensusQuantile} rather than \code{consensusQuantile},
be reproduced? Use this only to reproduce old calculations.}

  \item{networkOptions}{An object of class \code{\link{NetworkOptions}} giving the network construction
options to be used in the calculation of the similarity.}

  \item{returnDissim}{Logical: if \code{TRUE}, dissimarity, rather than similarity, will be returned.}

  \item{consensusTree}{A list of class \code{\link{ConsensusTree}} specifying the consensus calculation.
Note that calibration options within the
consensus specifications are ignored: since the consensus is calulated from entries representing a single
value, calibration would not make sense.}

  \item{\dots}{
Other arguments for compatibility; currently unused.
}
}
\details{
These functions calculate the similarity or dissimilarity of two groups of genes (variables) in \code{expr} or
\code{multiExpr} using correlations of the first singular vectors ("eigengenes"). For a single data set
(\code{branchEigengeneDissim} and \code{branchEigengeneSimilarity}), the similarity is the correlation, and
dissimilarity 1-correlation of the first signular vectors. 

Functions \code{mtd.branchEigengeneDissim} and 
\code{hierarchicalBranchEigengeneDissim} calculate consensus eigengene dissimilarity.
Function \code{mtd.branchEigengeneDissim} calculates a simple ("flat") consensus of branch eigengene
similarities across the given data set, at the given consensus quantile. 
Function \code{hierarchicalBranchEigengeneDissim} can calculate a hierarchical consensus in which consensus
calculations are hierarchically nested. 
}
\value{
A single number, the dissimilarity for \code{branchEigengeneDissim}, \code{mtd.branchEigengeneDissim}, and
\code{hierarchicalBranchEigengeneDissim}. 

\code{branchEigengeneSimilarity} returns similarity or dissimilarity, depending on imput.
}
\author{
Peter Langfelder
}
\seealso{\code{\link{hierarchicalConsensusCalculation}}}
\keyword{ misc }% __ONLY ONE__ keyword per line
