\name{consensusProjectiveKMeans}
\alias{consensusProjectiveKMeans}
\title{ Consensus projective K-means (pre-)clustering of expression data }
\description{
  Implementation of a consensus variant of K-means clustering for expression data across multiple data sets.
}
\usage{
consensusProjectiveKMeans(
  multiExpr, 
  preferredSize = 5000, 
  nCenters = NULL, 
  sizePenaltyPower = 4,
  networkType = "unsigned", 
  randomSeed = 54321,
  checkData = TRUE,
  imputeMissing = TRUE,
  useMean = (length(multiExpr) > 3),
  maxIterations = 1000, 
  verbose = 0, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiExpr}{  expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{preferredSize}{ preferred maximum size of clusters. }

  \item{nCenters}{ number of initial clusters. Empirical evidence suggests that more centers will give a 
better preclustering; the default is \code{as.integer(min(nGenes/20, 100*nGenes/preferredSize))} 
and is an attempt to arrive at a reasonable number given the resources available. }

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
start? If the missing data are not imputed, they will be replaced by 0 which can be problematic.}


  \item{useMean}{ logical: should mean distance across sets be used instead of maximum? See details. }

  \item{maxIterations}{ maximum iterations to be attempted. }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{

The principal aim of this function within WGCNA is to pre-cluster a large number of genes into smaller
blocks that can be handled using standard WGCNA techniques. 

This function implements a variant of K-means clustering that is suitable for co-expression analysis.
Cluster centers are defined by the first principal component, and distances by correlation. Consensus
distance across several sets is defined as the maximum of the corresponding distances in individual
sets; however, if \code{useMean} is set, the mean distance will be used instead of the maximum. 
The distance between a gene and a center of a cluster is multiplied by a factor of
\eqn{max(clusterSize/preferredSize, 1)^{sizePenaltyPower}}{\code{max(clusterSize/preferredSize,
1)^sizePenaltyPower}}, thus penalizing clusters whose size exceeds
\code{preferredSize}. The function starts with randomly generated cluster assignment (hence the need to
set the random seed for repeatability) and executes interations of calculating new centers and
reassigning genes to nearest (in the consensus sense) center until the clustering becomes stable. 
Before returning, nearby
clusters are iteratively combined if their combined size is below \code{preferredSize}.

Consensus distance defined as maximum of distances in all sets is consistent with the approach taken in
\code{\link{blockwiseConsensusModules}}, but the procedure may not converge. Hence it is advisable to use
the mean as consensus in cases where there are multiple data sets (4 or more, say) and/or if the input
data sets are very different. 

The standard principal component calculation via the function \code{svd} fails from time to time
(likely a convergence problem of the underlying lapack functions). Such errors are trapped and the
principal component is approximated by a weighted average of expression profiles in the cluster. If
\code{verbose} is set above 2, an informational message is printed whenever this approximation is used.
}
 
\value{
 A list with the following components:
 \item{clusters}{ a numerical vector with one component per input gene, giving the cluster number in
which the gene is assigned. }

 \item{centers}{ a vector of lists, one list per set. Each list contains a component \code{data} that
contains a matrix whose columns are the cluster centers in the corresponding set. }

  \item{unmergedClusters}{  a numerical vector with one component per input gene, giving the cluster
number in which the gene was assigned before the final merging step. }

  \item{unmergedCenters}{  a vector of lists, one list per set. Each list contains a component
\code{data} that contains a matrix whose columns are the cluster centers before merging in the
corresponding set. } 

}
\author{ Peter Langfelder }
\seealso{ \code{\link{projectiveKMeans}} }
\keyword{ cluster }
