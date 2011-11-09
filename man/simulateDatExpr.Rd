\name{simulateDatExpr}
\alias{simulateDatExpr}
\title{ Simulation of expression data}
\description{
  Simulation of expression data with a customizable modular structure and several different types of
noise. 
}
\usage{
simulateDatExpr(
  eigengenes, 
  nGenes, 
  modProportions, 
  minCor = 0.3, 
  maxCor = 1, 
  corPower = 1, 
  signed = FALSE, 
  propNegativeCor = 0.3, 
  geneMeans = NULL,
  backgroundNoise = 0.1, 
  leaveOut = NULL, 
  nSubmoduleLayers = 0, 
  nScatteredModuleLayers = 0, 
  averageNGenesInSubmodule = 10, 
  averageExprInSubmodule = 0.2, 
  submoduleSpacing = 2, 
  verbose = 1, indent = 0)
}

\arguments{
  \item{eigengenes}{ a data frame containing the seed eigengenes for the simulated modules. Rows
correspond to samples and columns to modules. }

  \item{nGenes}{ total number of genes to be simulated. }

  \item{modProportions}{ a numeric vector with length equal the number of eigengenes in \code{eigengenes}
plus one, containing fractions of the total number of genes to be put into each of the modules and into
the "grey module", which means genes not related to any of the modules. See details. }

  \item{minCor}{ minimum correlation of module genes with the corresponding eigengene. See details. }

  \item{maxCor}{ maximum correlation of module genes with the corresponding eigengene. See details. }

  \item{corPower}{ controls the dropoff of gene-eigengene correlation. See details. }

  \item{signed}{ logical: should the genes be simulated as belonging to a signed network? If \code{TRUE},
all genes will be simulated to have positive correlation with the eigengene. If \code{FALSE}, a
proportion given by \code{propNegativeCor} will be simulated with negative correlations of the same
absolute values. }

  \item{propNegativeCor}{ proportion of genes to be simulated with negative gene-eigengene correlations.
Only effective if \code{signed} is \code{FALSE}. }

  \item{geneMeans}{ optional vector of length \code{nGenes} giving desired mean expression for each gene. If
not given, the returned expression profiles will have mean zero. }
 
  \item{backgroundNoise}{ amount of background noise to be added to the simulated expression data. }

  \item{leaveOut}{ optional specification of modules that should be left out of the simulation, that is
their genes will be simulated as unrelated ("grey"). This can
be useful when simulating several sets, in some which a module is present while in others it is absent. }

  \item{nSubmoduleLayers}{ number of layers of ordered submodules to be added. See details. }

  \item{nScatteredModuleLayers}{ number of layers of scattered submodules to be added. See details. }

  \item{averageNGenesInSubmodule}{ average number of genes in a submodule. See details. }

  \item{averageExprInSubmodule}{ average strength of submodule expression vectors. }
 
  \item{submoduleSpacing}{ a number giving submodule spacing: this multiple of the submodule size will
lie between the submodule and the next one.  }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{

  Given \code{eigengenes} can be unrelated or they can exhibit non-trivial correlations. Each module is
simulated separately from others. The expression profiles are chosen such that their
correlations with the eigengene run from just below \code{maxCor} to \code{minCor} (hence minCor must be
between 0 and 1, not including the bounds). The parameter \code{corPower} can be chosen to control the
behaviour of the simulated correlation with the gene index; values higher than 1 will result in the
correlation approaching \code{minCor} faster and lower than 1 slower. 

  Numbers of genes in each module are specified (as fractions of the total number of genes \code{nGenes})
by \code{modProportions}. The last entry in \code{modProportions} corresponds to the genes that will be
simulated as unrelated to anything else ("grey" genes). The proportion must add up to 1 or less. If the
sum is less than one, the remaining genes will be partitioned into groups and simulated to be "close" to
the proper modules, that is with small but non-zero correlations (between \code{minCor} and 0)
with the module eigengene. 

If \code{signed} is set \code{FALSE}, the correlation for
some of the module genes is chosen negative (but the absolute values remain the same as they would be for
positively correlated genes). To ensure consistency for simulations of multiple sets, the indices of the
negatively correlated genes are fixed and distributed evenly. 

In addition to the primary module structure, a secondary structure can be optionally simulated. Modules
in the secondary structure have sizes chosen from an exponential distribution with mean equal
\code{averageNGenesInSubmodule}. Expression vectors simulated in the secondary structure are simulated
with expected standard deviation chosen from an exponential distribution with mean equal
\code{averageExprInSubmodule}; the higher this coefficient, the
more pronounced will the submodules be in the main modules. The secondary structure can be simulated in
several layers; their number is given by \code{SubmoduleLayers}. Genes in these submodules are ordered in
the same order as in the main modules. 

In addition to the ordered submodule structure, a scattered submodule structure can be simulated as well.
This structure can be viewed as noise that tends to correlate random groups of genes. The size and effect 
parameters are the same as for the ordered submodules, and the number of layers added is controlled by
\code{nScatteredModuleLayers}. 
}
\value{
  A list with the following components:
  \item{datExpr}{ simulated expression data in a data frame whose columns correspond genes and rows to
samples. }

  \item{setLabels}{ simulated module assignment. Module labels are numeric, starting from 1. Genes
simulated to be outside of proper modules have label 0. 
Modules that are left out (specified in \code{leaveOut})
are indicated as 0 here. }

  \item{allLabels}{ simulated module assignment. Genes that belong to leftout modules (specified in
\code{leaveOut}) are indicated by their would-be  assignment here. }

  \item{labelOrder}{ a vector specifying the order in which labels correspond to the given eigengenes,
that is \code{labelOrder[1]} is the label assigned to module whose seed is \code{eigengenes[, 1]} etc.  }
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

\code{\link{simulateModule}} for simulations of individual modules;

\code{\link{simulateDatExpr5Modules}} for a simplified interface to expression simulations;

\code{\link{simulateMultiExpr}} for a simulation of several related data sets. 

}
\keyword{ misc }% __ONLY ONE__ keyword per line
