\name{nearestNeighborConnectivity}
\alias{nearestNeighborConnectivity}
\title{ Connectivity to a constant number of nearest neighbors }
\description{
  Given expression data and basic network parameters, the function calculates connectivity of each gene
to a given number of nearest neighbors.
}
\usage{
nearestNeighborConnectivity(datExpr, 
         nNeighbors = 50, power = 6, type = "unsigned", 
         corFnc = "cor", corOptions = "use = 'p'", 
         blockSize = 1000,
         sampleLinks = NULL, nLinks = 5000, setSeed = 38457,
         verbose = 1, indent = 0)
}
\arguments{
  \item{datExpr}{ a data frame containing expression data, with rows corresponding to samples and columns
to genes. Missing values are allowed and will be ignored. }

  \item{nNeighbors}{ number of nearest neighbors to use. }

  \item{power}{ soft thresholding power for network construction. Should be a number greater than 1. }

  \item{type}{ a character string encoding network type. Recognized values are (unique abbreviations of) 
\code{"unsigned"}, \code{"signed"}, and \code{"signed hybrid"}. }

  \item{corFnc}{ character string containing the name of the function to calculate correlation. Suggested
functions include \code{"cor"} and \code{"bicor"}. }

  \item{corOptions}{ further argument to the correlation function. }

  \item{blockSize}{ correlation calculations will be split into square blocks of this size, to prevent
running out of memory for large gene sets.  }

  \item{sampleLinks}{ logical: should network connections be sampled (\code{TRUE}) or should all
connections be used systematically (\code{FALSE})? }

  \item{nLinks}{ number of links to be sampled. Should be set such that \code{nLinks * nNeighbors} be
several times larger than the number of genes. }

  \item{setSeed}{ seed to be used for sampling, for repeatability. If a seed already exists, it is saved
before the sampling starts and restored upon exit. }

  \item{verbose}{ integer controlling the level of verbosity. 0 means silent.}

  \item{indent}{ integer controlling indentation of output. Each unit above 0 adds two spaces. }
}
\details{
  Connectivity of gene \code{i} is the sum of adjacency strengths between gene \code{i} 
and other genes; in
this case we take the \code{nNeighbors} nodes with the highest connection strength to gene \code{i}. The
adjacency strengths are calculated by correlating the given expression data using the function supplied
in \code{corFNC} and transforming them into adjacency according to the given network \code{type} and
\code{power}. 
}
\value{
  A vector with one component for each gene containing the nearest neighbor connectivity. 
}
\author{ Peter Langfelder }
\seealso{ \code{\link{adjacency}}, \code{\link{softConnectivity}} }
\keyword{ misc }% __ONLY ONE__ keyword per line
