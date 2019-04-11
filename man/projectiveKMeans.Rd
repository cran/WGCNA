\name{projectiveKMeans}
\alias{projectiveKMeans}
\title{ Projective K-means (pre-)clustering of expression data }
\description{
  Implementation of a variant of K-means clustering for expression data. 
}
\usage{
projectiveKMeans(
  datExpr, 
  preferredSize = 5000, 
  nCenters = as.integer(min(ncol(datExpr)/20, preferredSize^2/ncol(datExpr))),
  sizePenaltyPower = 4, 
  networkType = "unsigned",  
  randomSeed = 54321,
  checkData = TRUE,
  imputeMissing = TRUE,
  maxIterations = 1000, 
  verbose = 0, indent = 0)
}
\arguments{

  \item{datExpr}{  expression data. A data frame in which columns are genes and rows ar samples. NAs are
allowed, but not too many. }

  \item{preferredSize}{ preferred maximum size of clusters. }

  \item{nCenters}{ number of initial clusters. Empirical evidence suggests that more centers will give a
better preclustering; the default is an attempt to arrive at a reasonable number. }

  \item{sizePenaltyPower}{ parameter specifying how severe is the penalty for clusters that exceed
\code{preferredSize}. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{randomSeed}{ integer to be used as seed for the random number generator before the function
starts. If a current seed exists, it is saved and restored upon exit. }

  \item{checkData}{ logical: should data be checked for genes with zero variance and 
genes and samples with excessive numbers of missing samples? Bad samples are ignored; returned cluster
assignment for bad genes will be \code{NA}. }

  \item{imputeMissing}{ logical: should missing values in \code{datExpr} be imputed before the calculations
start? The early imputation makes the code run faster but may produce slightly different results if re-running
older calculations.}

  \item{maxIterations}{ maximum iterations to be attempted. }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{

The principal aim of this function within WGCNA is to pre-cluster a large number of genes into smaller blocks
that can be handled using standard WGCNA techniques.

This function implements a variant of K-means clustering that is suitable for co-expression analysis.
Cluster centers are defined by the first principal component, and distances by correlation (more
precisely, 1-correlation). The distance between a gene and a cluster is multiplied by a factor of
\eqn{max(clusterSize/preferredSize, 1)^{sizePenaltyPower}}{\code{max(clusterSize/preferredSize,
1)^sizePenaltyPower}}, thus penalizing clusters whose size exceeds
\code{preferredSize}. The function starts with randomly generated cluster assignment (hence the need to
set the random seed for repeatability) and executes interations of calculating new centers and
reassigning genes to nearest center until the clustering becomes stable. Before returning, nearby
clusters are iteratively combined if their combined size is below \code{preferredSize}.

The standard principal component calculation via the function \code{svd} fails from time to time
(likely a convergence problem of the underlying lapack functions). Such errors are trapped and the
principal component is approximated by a weighted average of expression profiles in the cluster. If
\code{verbose} is set above 2, an informational message is printed whenever this approximation is used.

}
\value{
  A list with the following components:
  \item{clusters}{A numerical vector with one component per input gene, giving the cluster number in
which the gene is assigned. }
  \item{centers}{Cluster centers, that is their first principal components. }
}
\author{ Peter Langfelder }
\seealso{ \code{\link{sizeRestrictedClusterMerge}} which implements the last step of merging smaller clusters.}
\keyword{ cluster }
