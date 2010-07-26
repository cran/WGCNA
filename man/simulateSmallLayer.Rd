\name{simulateSmallLayer}
\alias{simulateSmallLayer}
\title{ Simulate small modules }
\description{
  This function simulates a set of small modules. The primary purpose is to add a submodule structure to
the main module structure simulated by \code{\link{simulateDatExpr}}. 
}
\usage{
simulateSmallLayer(
  order, 
  nSamples, 
  minCor = 0.3, maxCor = 0.5, corPower = 1, 
  averageModuleSize, 
  averageExpr, 
  moduleSpacing, 
  verbose = 4, indent = 0)
}
\arguments{
  \item{order}{ a vector giving the simulation order for vectors. See details. }

  \item{nSamples}{ integer giving the number of samples to be simulated. }

  \item{minCor}{ a multiple of \code{maxCor} (see below) giving the minimum correlation of module genes
with the corresponding eigengene. See details. } 

  \item{maxCor}{ maximum correlation of module genes with the corresponding eigengene. See details. }

  \item{corPower}{ controls the dropoff of gene-eigengene correlation. See details. }

  \item{averageModuleSize}{  average number of genes in a module. See details. }

  \item{averageExpr}{ average strength of module expression vectors. }

  \item{moduleSpacing}{ a number giving module spacing: this multiple of the module size will
lie between the module and the next one.  }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
  Module eigenvectors are chosen randomly and independently. Module sizes are chosen randomly from an
exponential distribution with mean equal \code{averageModuleSize}. Two thirds of genes in 
each module are simulated
as proper module genes and one third as near-module genes (see \code{\link{simulateModule}} for details). 
Between each successive pairs of
modules a number of genes given by \code{moduleSpacing} will be left unsimulated (zero expression).
Module expression, that is the
expected standard deviation of the module expression vectors, is chosen randomly from an exponential
distribution with mean equal \code{averageExpr}. The expression profiles are chosen such that their
correlations with the eigengene run from just below \code{maxCor} to \code{minCor * maxCor} 
(hence minCor must be
between 0 and 1, not including the bounds). The parameter \code{corPower} can be chosen to control the
behaviour of the simulated correlation with the gene index; values higher than 1 will result in the
correlation approaching \code{minCor * maxCor} faster and lower than 1 slower. 

The simulated genes will be returned in the order given in \code{order}. 
}
\value{
 A matrix of simulated gene expressions, with dimension \code{(nSamples, length(order))}. 
}
\author{ Peter Langfelder }

\seealso{ 
 
  \code{\link{simulateModule}} for simulation of individual modules;

  \code{\link{simulateDatExpr}} for the main gene expression simulation function.
}

\keyword{ misc }% __ONLY ONE__ keyword per line
