\name{preservationNetworkConnectivity}
\alias{preservationNetworkConnectivity}
\title{ Network preservation calculations }
\description{
  This function calculates several measures of gene network preservation. Given gene expression data in
several individual data sets, it calculates the individual adjacency matrices, forms the preservation
network and finally forms several summary measures of adjacency preservation for each node (gene) in the
network.
}
\usage{
preservationNetworkConnectivity(
   multiExpr,
   useSets = NULL, useGenes = NULL,
   corFnc = "cor", corOptions = "use='p'",
   networkType = "unsigned",
   power = 6,
   sampleLinks = NULL, nLinks = 5000,
   blockSize = 1000,
   setSeed = 12345,
   weightPower = 2,
   verbose = 2, indent = 0)

}
\arguments{
  \item{multiExpr}{ expression data in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{useSets}{ optional specification of sets to be used for the preservation calculation. Defaults to
using all sets. }

  \item{useGenes}{ optional specification of genes to be used for the preservation calculation. Defaults
to all genes. }

 \item{corFnc}{ character string containing the name of the function to calculate correlation. Suggested
functions include \code{"cor"} and \code{"bicor"}. }

  \item{corOptions}{ further argument to the correlation function. }

  \item{networkType}{ a character string encoding network type. Recognized values are (unique
abbreviations of) \code{"unsigned"}, \code{"signed"}, and \code{"signed hybrid"}. }

  \item{power}{ soft thresholding power for network construction. Should be a number greater than 1. }

  \item{sampleLinks}{ logical: should network connections be sampled (\code{TRUE}) or should all
connections be used systematically (\code{FALSE})? }

  \item{nLinks}{ number of links to be sampled. Should be set such that \code{nLinks * nNeighbors} be
several times larger than the number of genes. }

  \item{blockSize}{ correlation calculations will be split into square blocks of this size, to prevent
running out of memory for large gene sets.  }

  \item{setSeed}{  seed to be used for sampling, for repeatability. If a seed already exists, it is saved
before the sampling starts and restored upon exit. }

  \item{weightPower}{ power with which higher adjacencies will be weighted in weighted means }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
The preservation network is formed from adjacencies of compared sets. For 'complete' preservations, all
given sets are compared at once; for 'pairwise' preservations, the sets are compared in pairs. Unweighted
preservations are simple mean preservations for each node; their weighted counterparts are weighted
averages in which a preservation of adjacencies \eqn{A^{(1)}_{ij}}{A[i,j; 1]} and 
\eqn{A^{(2)}_{ij}}{A[i,j; 2]} of nodes \eqn{i,j} between sets 1 and 2
is weighted by 
\eqn{[ (A^{(1)}_{ij} + A^{(2)}_{ij} )/2]^weightPower}{ ( (A[i,j; 1]+A[i,j; 2])/2)^weightPower}.
The hyperbolic preservation is based on 
\eqn{tanh[( max - min)/(max+min)^2]}{tanh[( max - min)/(max+min)^2]}, where \eqn{max}{max} and
\eqn{min}{min} are the componentwise maximum and minimum of the compared adjacencies, respectively.
  
}
\value{
  A list with the following components:

  \item{pairwise}{ a matrix with rows corresponding to genes and columns to unique pairs of given sets,
giving the pairwise preservation of the adjacencies connecting the gene to all other genes.}

  \item{complete}{ a vector with one entry for each input gene containing the complete mean preservation
of the adjacencies connecting the gene to all other genes.}

  \item{pairwiseWeighted}{ a matrix with rows corresponding to genes and columns to unique pairs of given
sets, giving the pairwise weighted preservation of the adjacencies connecting the gene to all other genes.}

  \item{completeWeighted}{ a vector with one entry for each input gene containing the complete weighted 
mean preservation of the adjacencies connecting the gene to all other genes.}

  \item{pairwiseHyperbolic}{ a matrix with rows corresponding to genes and columns to unique 
pairs of given sets, giving the pairwise hyperbolic preservation of the adjacencies connecting the gene
to all other genes.}

  \item{completeHyperbolic}{ a vector with one entry for each input gene containing the complete mean
hyperbolic preservation of the adjacencies connecting the gene to all other genes.}

  \item{pairwiseWeightedHyperbolic}{ a matrix with rows corresponding to genes and columns to unique
pairs of given sets, giving the pairwise weighted hyperbolic preservation of the adjacencies connecting
the gene to all other genes.} 

  \item{completeWeightedHyperbolic}{ a vector with one entry for each input gene containing the complete
weighted hyperbolic  mean preservation of the adjacencies connecting the gene to all other genes.}
}
\references{ Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between
co-expression modules. BMC Systems Biology 2007, 1:54 }

\author{ Peter Langfelder }

\seealso{ 
  \code{\link{adjacency}} for calculation of adjacency; 
}
\keyword{misc}
