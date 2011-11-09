\name{simulateMultiExpr}
\alias{simulateMultiExpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate multi-set expression data}

\description{
  Simulation of expression data in several sets with relate module structure.
}
\usage{
simulateMultiExpr(eigengenes, 
                  nGenes, 
                  modProportions, 
                  minCor = 0.5, maxCor = 1, 
                  corPower = 1, 
                  backgroundNoise = 0.1, 
                  leaveOut = NULL, 
                  signed = FALSE, 
                  propNegativeCor = 0.3, 
                  geneMeans = NULL,
                  nSubmoduleLayers = 0, 
                  nScatteredModuleLayers = 0, 
                  averageNGenesInSubmodule = 10, 
                  averageExprInSubmodule = 0.2, 
                  submoduleSpacing = 2, 
                  verbose = 1, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eigengenes}{  the seed eigengenes for the simulated modules in a multi-set format. A list with one
component per set. Each component is again a list that must contain a component \code{data}. This is a data
frame of seed eigengenes for the corresponding data set. Columns correspond to modules, rows to samples.
Number of samples in the simulated data is determined from the number of samples of the eigengenes. }

  \item{nGenes}{ integer specifyin the number of simulated genes. }

  \item{modProportions}{  a numeric vector with length equal the number of eigengenes in \code{eigengenes}
plus one, containing fractions of the total number of genes to be put into each of the modules and into
the "grey module", which means genes not related to any of the modules. See details. }

  \item{minCor}{ minimum correlation of module genes with the corresponding eigengene. See details. }

  \item{maxCor}{ maximum correlation of module genes with the corresponding eigengene. See details. }

  \item{corPower}{ controls the dropoff of gene-eigengene correlation. See details. }

  \item{backgroundNoise}{ amount of background noise to be added to the simulated expression data. }

    \item{leaveOut}{ optional specification of modules that should be left out of the simulation, that is
their genes will be simulated as unrelated ("grey"). A logical matrix in which columns correspond to sets
and rows to modules. Wherever \code{TRUE}, the corresponding module in the corresponding data set will not
be simulated, that is its genes will be simulated independently of the eigengene. }

   \item{signed}{ logical: should the genes be simulated as belonging to a signed network? If \code{TRUE},
all genes will be simulated to have positive correlation with the eigengene. If \code{FALSE}, a
proportion given by \code{propNegativeCor} will be simulated with negative correlations of the same
absolute values. }
 
  \item{propNegativeCor}{ proportion of genes to be simulated with negative gene-eigengene correlations.
Only effective if \code{signed} is \code{FALSE}. }

 \item{geneMeans}{ optional vector of length \code{nGenes} giving desired mean expression for each gene. If
not given, the returned expression profiles will have mean zero. }

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
  For details of simulation of individual data sets and the meaning of individual set simulation arguments, 
see \code{\link{simulateDatExpr}}. This function
simulates several data sets at a time and puts the result in a multi-set format. The number of genes is the
same for all data sets. Module memberships are also the same, but modules can optionally be ``dissolved'',
that is their genes will be simulated as unassigned. Such ``dissolved'', or left out, modules can be
specified in the matrix \code{leaveOut}.
}
\value{
  A list with the following components:

  \item{multiExpr }{simulated expression data in multi-set format analogous to that of the input
\code{eigengenes}.  A list with one
component per set. Each component is again a list that must contains a component \code{data}. This is a data
frame of expression data for the corresponding data set. Columns correspond to genes, rows to samples.}

  \item{setLabels}{a matrix of dimensions (number of genes) times (number of sets) that contains module
labels for each genes in each simulated data set. }

 \item{allLabels}{a matrix of dimensions (number of genes) times (number of sets) that contains the module
labels that would be simulated if no module were left out using \code{leaveOut}. This means that all columns
of the matrix are equal; the columns are repeated for convenience so \code{allLabels} has the same
dimensions as \code{setLabels}. }

 \item{labelOrder}{a matrix of dimensions (number of modules) times (number of sets) that contains the 
order in which module labels were assigned to genes in each set. The first label is assigned to genes
1...(module size of module labeled by first label), the second label to the following batch of genes etc.} 

}

\references{ 
 A short description of the simulation method can also be found in the Supplementary Material
to the article

Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between co-expression
modules. BMC Systems Biology 2007, 1:54.

The material is posted at
http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/EigengeneNetwork/SupplementSimulations.pdf.

}
\author{ Peter Langfelder}
\seealso{ 

\code{\link{simulateEigengeneNetwork}} for a simulation of eigengenes with a given causal structure;

\code{\link{simulateDatExpr}} for simulation of individual data sets; 

\code{\link{simulateDatExpr5Modules}} for a simple simulation of a data set consisting of 5 modules;

\code{\link{simulateModule}} for simulations of individual modules;

}
\keyword{ misc }% __ONLY ONE__ keyword per line
