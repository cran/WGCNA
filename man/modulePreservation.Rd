\name{modulePreservation}
\alias{modulePreservation}
\title{ Calculation of module preservation statistics }
\description{
  Calculations of module preservation statistics between independent data sets.
}
\usage{
modulePreservation(
   multiData,
   multiColor,
   dataIsExpr = TRUE,
   networkType = "unsigned", 
   corFnc = "cor",
   corOptions = "use = 'p'",
   referenceNetworks = 1, 
   testNetworks = NULL,
   nPermutations = 100, 
   includekMEallInSummary = FALSE,
   restrictSummaryForGeneralNetworks = TRUE,
   calculateQvalue = FALSE,
   randomSeed = 12345, 
   maxGoldModuleSize = 1000, 
   maxModuleSize = 1000, 
   quickCor = 1, 
   ccTupletSize = 2, 
   calculateCor.kIMall = FALSE,
   calculateClusterCoeff = FALSE,
   useInterpolation = FALSE, 
   checkData = TRUE, 
   greyName = NULL, 
   savePermutedStatistics = TRUE, 
   loadPermutedStatistics = FALSE, 
   permutedStatisticsFile = if (useInterpolation) "permutedStats-intrModules.RData" 
                                   else "permutedStats-actualModules.RData", 
   plotInterpolation = TRUE, 
   interpolationPlotFile = "modulePreservationInterpolationPlots.pdf", 
   discardInvalidOutput = TRUE,
   parallelCalculation = FALSE,
   verbose = 1, indent = 0)
}
\arguments{
  \item{multiData}{  expression data or adjacency data 
in the multi-set format (see \code{\link{checkSets}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression or adjacency
data.
If expression data are used,
rows correspond to samples and columns to genes or probes. In case of adjacencies, each \code{data} matrix
should be a symmetric matrix ith entries between 0 and 1 and unit diagonal. 
Each component of the outermost list should be
named. }

  \item{multiColor}{ a list in which every component is a vector giving the module labels of genes in
\code{multiExpr}. The components must be named using the same names that are used in \code{multiExpr}; these
names are used top match labels to expression data sets. See details. }

  \item{dataIsExpr}{ logical: if \code{TRUE}, \code{multiData} will be interpreted as expression data; if
\code{FALSE}, \code{multiData} will be interpreted as adjacencies. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{corFnc}{ character string specifying the function to be used to calculate co-expression
similarity. Defaults to Pearson correlation. Another useful choice is \code{\link{bicor}}.
More generally, any function returning values between -1 and 1 can be used. }

  \item{corOptions}{ character string specifying additional arguments to be passed to the function given
by \code{corFnc}. Use \code{"use = 'p', method = 'spearman'"} to obtain Spearman correlation.   }

  \item{referenceNetworks}{ a vector giving the indices of expression data to be used as reference networks.
Reference networks must have their module labels given in \code{multiColor}. }

  \item{testNetworks}{a list with one component per each entry in \code{referenceNetworks} above, giving
the test networks in which to evaluate module preservation for the corresponding reference network. If not
given, preservation will be evaluated in all networks (except each reference network). If
\code{referenceNetworks} is of length 1, \code{testNetworks} can also be a vector (instead of a list
containing the single vector).} 

  \item{nPermutations}{ specifies the number of permutations that will be calculated in the permutation test. }

  \item{includekMEallInSummary}{ logical: should cor.kMEall be included in the calculated summary statistics?
Because kMEall takes into account all genes in the network, this statistic measures preservation of the full
network with respect to the eigengene of the module. This may be undesirable, hence the default is
\code{FALSE}.}

  \item{restrictSummaryForGeneralNetworks}{ logical: should the summary statistics for general (not
correlation) networks be restricted (density to meanAdj, connectivity to cor.kIM and cor.Adj)? The default
\code{TRUE} corresponds to published work. }

  \item{calculateQvalue}{ logical: should q-values (local FDR estimates) be calculated? Package qvalue must
be installed for this calculation. Note that q-values may not be meaningful when the number of modules is
small and/or most modules are preserved. }

  \item{randomSeed}{ seed for the random number generator. If \code{NULL}, the seed will not be set. If
non-\code{NULL} and the random generator has been initialized prior to the function call, the latter's state
is saved and restored upon exit} 

  \item{maxGoldModuleSize}{ maximum size of the "gold" module, i.e., the random sample of all network genes. }

  \item{maxModuleSize}{ maximum module size used for calculations. Modules larger than \code{maxModuleSize}
will be reduced by randomly sampling \code{maxModuleSize} genes. }

  \item{quickCor}{ number between 0 and 1 specifying the handling of missing data in calculation of
correlation. Zero means exact but potentially slower calculations; one means potentially faster
calculations, but with potentially inaccurate results if the proportion of missing data is large. See
\code{\link{cor}} for more details. }

  \item{ccTupletSize}{ tuplet size for co-clustering calculations. }

  \item{calculateCor.kIMall}{ logical: should cor.kMEall be calculated? This option is only valid for
adjacency input. If \code{FALSE}, cor.kIMall will not be calculated, potentially saving significant amount
of time if the input adjacencies are large and contain many modules. }

  \item{calculateClusterCoeff}{ logical: should statistics based on the clustering coefficient be
calculated? While these statistics may be interesting, the calculations are also computationally expensive.}

  \item{checkData}{ logical: should data be checked for excessive number of missing entries? See
\code{\link{goodSamplesGenesMS}} for details. }

  \item{greyName}{ label used for unassigned genes. Traditionally such genes are labeled by grey color or
numeric label 0. These values are the default when \code{multiColor} contains character or numeric vectors,
respectively. }

  \item{savePermutedStatistics}{ logical: should calculated permutation statistics be saved? Saved
statistics may be re-used if the calculation needs to be repeated.}

  \item{permutedStatisticsFile}{ file name to save the permutation statistics into. }

  \item{loadPermutedStatistics}{ logical: should permutation statistics be loaded? If a previously executed
calculation needs to be repeated, loading permutation study results can cut the calculation time many-fold. }

  \item{useInterpolation}{ logical: should permutation statistics be calculated by interpolating an artificial
set of evenly spaced modules? This option may potentially speed up the calculations, but it restricts
calculations to density measures.  }

  \item{plotInterpolation}{ logical: should interpolation plots be saved? If interpolation is used (see
\code{useInterpolation} above), the function can optionally generate diagnostic plots that can be used to
assess whether the interpolation makes sense. }

  \item{interpolationPlotFile}{ file name to save the interpolation plots into. }

  \item{discardInvalidOutput}{logical: should output columns containing no valid data be discarded? This
option may be useful when input \code{dataIsExpr} is \code{FALSE} and some of the output statistics cannot
be calculated. This option causes such statistics to be dropped from output.}

  \item{parallelCalculation}{logical: should calculations be done in parallel? Note that parallel
calculations are turned off by default and will lead to somewhat DIFFERENT results than serial calculations
because the random seed is set differently. For the calculation to actually run in parallel mode, a call to
\code{\link{enableWGCNAThreads}} must be made before this function is called.} 

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
  This function calculates module preservation statistics pair-wise between given reference sets and all
other sets in \code{multiExpr}. Reference sets must have their corresponding module assignment specified in
\code{multiColor}; module assignment is optional for test sets. Individual expression sets and their module
labels are matched using \code{names} of the corresponding components in \code{multiExpr} and
\code{multiColor}. 

For each reference-test pair, the function calculates module preservation statistics that
measure how well the modules of the reference set are preserved in the test set. 
If the \code{multiColor} also contains module assignment for the test set, the calculated statistics also
include cross-tabulation statistics that make use of the test module assignment. 

For each reference-test pair, the function only uses genes (columns of the \code{data} component of each
component of \code{multiExpr}) that are in common between the reference and test set. Columns are matched by
column names, so column names must be valid.  

In addition to preservation statistics, the function also calculates several statistics of module quality,
that is measures of how well-defined modules are in the reference set. The quality statistics are calculated
with respect to genes in common with with a test set; thus the function calculates a set of quality
statistics for each reference-test pair. This may be somewhat counter-intuitive, but it allows a direct
comparison of corresponding quality and preservation statistics.

The calculated p-values are determined from the Z scores of individual measures under assumption of
normality. No p-value is calculated for the Zsummary measures. Bonferoni correction to the number of tested
modules. Because the p-values for strongly preserved modules are often extremely low, the function reports
natural logarithms (base e) of the p-values. However, q-values are reported untransformed since they are
calculated that way in package qvalue.

Missing data are removed (but see \code{quickCor} above). 

}
\value{
  The function returns a nested list of preservation statistics. At the top level, the list components are:
  \item{quality}{observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally)
q-values of quality statistics. All logarithms are in base 10.}
  \item{preservation }{observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and
(optionally) q-values of density and connectivity preservation statistics. All logarithms are in base 10.}
  \item{accuracy}{observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and
(optionally) q-values of cross-tabulation statistics. All logarithms are in base 10.}
  \item{referenceSeparability}{observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally)
q-values of module separability in the reference network. All logarithms are in base 10.}
  \item{testSeparability}{observed values, Z scores, p-values, Bonferoni-corrected p-values, and (optionally)
q-values of module separability in the test network. All logarithms are in base 10.}
  \item{permutationDetails}{results of individual permutations, useful for diagnostics}
  All of the above are lists. The lists \code{quality}, \code{preservation}, \code{referenceSeparability},
and \code{testSeparability} each contain 4 or 5 components: \code{observed} contains observed values, 
\code{Z} contains the corresponding Z scores, \code{log.p} contains base 10 logarithms of the p-values,
\code{log.pBonf}  contains base 10 logarithms of the Bonferoni corrected p-values, and optionally \code{q}
contains the associated q-values. The list \code{accuracy} contains \code{observed}, \code{Z}, \code{log.p},
\code{log.pBonf}, optionally \code{q}, 
and additional components \code{observedOverlapCounts} and \code{observedFisherPvalues} that contain the
observed matrices of overlap counts and Fisher test p-values. 

Each of the lists \code{observed}, \code{Z}, \code{log.p},
\code{log.pBonf}, optionally \code{q}, \code{observedOverlapCounts} and \code{observedFisherPvalues}
is structured as a 2-level list where the outer components correspond to reference sets and the inner
components to tests sets. As an example, \code{preservation$observed[[1]][[2]]} contains the density and
connectivity preservation statistics for the preservation of set 1 modules in set 2, that is set 1 is the
reference set and set 2 is the test set. \code{preservation$observed[[1]][[2]]} is a data frame in which
each row corresponds to a module in the reference network 1 plus one row for the unassigned objects, and
one row for a "module" that contains randomly sampled objects and that represents a whole-network average.
Each column corresponds to a statistic as indicated by the column name. 
}
\references{ Peter Langfelder, Rui Luo, Michael C. Oldham, and Steve Horvath, to appear }
\author{ Rui Luo and Peter Langfelder }
\note{ 

For large data sets, the permutation study may take a while (typically on the order of several hours). Use
 \code{verbose = 3} to get detailed progress report as the calculations advance.

}
\seealso{ Network construction and module detection functions in the WGCNA package such as
\code{\link{adjacency}}, \code{\link{blockwiseModules}}; rudimentary cleaning in
\code{\link{goodSamplesGenesMS}}; the WGCNA implementation of correlation in \code{\link{cor}}.

 }
\keyword{ misc }% __ONLY ONE__ keyword per line
