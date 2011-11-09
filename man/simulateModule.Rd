\name{simulateModule}
\alias{simulateModule}
\title{ Simulate a gene co-expression module}
\description{
  Simulation of a single gene co-expression module.
}
\usage{
simulateModule(
  ME, 
  nGenes, 
  nNearGenes = 0, 
  minCor = 0.3, maxCor = 1, corPower = 1, 
  signed = FALSE, propNegativeCor = 0.3, 
  geneMeans = NULL,
  verbose = 0, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ME}{ seed module eigengene. }

  \item{nGenes}{ number of genes in the module to be simulated. Must be non-zero. } 

  \item{nNearGenes}{ number of genes to be simulated with low correlation with the seed eigengene. }

  \item{minCor}{ minimum correlation of module genes with the eigengene. See details. }

  \item{maxCor}{ maximum correlation of module genes with the eigengene. See details. }

  \item{corPower}{  controls the dropoff of gene-eigengene correlation. See details.  }

  \item{signed}{ logical: should the genes be simulated as belonging to a signed network? If \code{TRUE},
all genes will be simulated to have positive correlation with the eigengene. If \code{FALSE}, a
proportion given by \code{propNegativeCor} will be simulated with negative correlations of the same
absolute values.  }

  \item{propNegativeCor}{ proportion of genes to be simulated with negative gene-eigengene correlations.
Only effective if \code{signed} is \code{FALSE}. }

  \item{geneMeans}{ optional vector of length \code{nGenes} giving desired mean expression for each gene. If
not given, the returned expression profiles will have mean zero. }

 \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
  Module genes are simulated around the eigengene by choosing them such that their (expected)
correlations with the seed eigengene decrease progressively from (just below) \code{maxCor} to \code{minCor}.
The genes are otherwise independent from one another. The variable \code{corPower} determines how fast
the correlation drops towards \code{minCor}. Higher powers lead to a faster frop-off; \code{corPower} must be
above zero but need not be integer.

If \code{signed} is \code{FALSE}, the genes are simulated so as to be part of an unsigned network module,
that is some genes will be simulated with a negative correlation with the seed eigengene (but of the same
absolute value that a positively correlated gene would be simulated with). The proportion of genes with
negative correlation is controlled by \code{propNegativeCor}. 

Optionally, the function can also simulate genes that are "near" the module, meaning they are
simulated with a  low but non-zero correlation with the seed eigengene. The correlations run between
\code{minCor} and zero. 

}
\value{  
  A matrix containing the expression data with rows corresponding to samples and columns to genes. 
}

\references{ A short description of the simulation method can also be found in the Supplementary Material
to the article

Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between co-expression
modules. BMC Systems Biology 2007, 1:54.

The material is posted at
http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/EigengeneNetwork/SupplementSimulations.pdf.
}

\author{ Peter Langfelder }
\seealso{ 

\code{\link{simulateEigengeneNetwork}} for a simulation of eigengenes with a given causal structure;

\code{\link{simulateDatExpr}} for simulations of whole datasets consisting of multiple modules;

\code{\link{simulateDatExpr5Modules}} for a simplified interface to expression simulations;

\code{\link{simulateMultiExpr}} for a simulation of several related data sets.

}
\keyword{ misc }% __ONLY ONE__ keyword per line
