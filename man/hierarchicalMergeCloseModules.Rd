\name{hierarchicalMergeCloseModules}
\alias{hierarchicalMergeCloseModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge close (similar) hierarchical consensus modules
}
\description{
Merges hierarchical consensus modules that are too close as measured by the correlation of their
eigengenes.
}
\usage{
hierarchicalMergeCloseModules(
  # input data
  multiExpr, 
  multiExpr.imputed = NULL,
  labels,

  # Optional starting eigengenes
  MEs = NULL,

  unassdColor = if (is.numeric(labels)) 0 else "grey",
  # If missing data are present, impute them?
  impute = TRUE,


  # Options for eigengene network construction
  networkOptions,

  # Options for constructing the consensus
  consensusTree,
  calibrateMESimilarities = FALSE,

  # Merging options
  cutHeight = 0.2,
  iterate = TRUE,

  # Output options
  relabel = FALSE,
  colorSeq = NULL,
  getNewMEs = TRUE,
  getNewUnassdME = TRUE,

  # Options controlling behaviour of the function
  trapErrors = FALSE,
  verbose = 1, indent = 0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{multiExpr}{ Expression data in the multi-set format (see \code{\link{multiData}}). A vector of
lists, one per set. Each set must contain a component \code{data} that contains the expression data, with
rows corresponding to samples and columns to genes or probes. }

  \item{multiExpr.imputed}{If \code{multiExpr} contain missing data, this argument can be used to supply the
expression data with missing data imputed. If not given, the \code{\link[impute]{impute.knn}} function will
be used to impute the missing data within each module (see \code{\link{imputeByModule}}.}

\item{labels}{
A vector (numeric, character or a factor) giving module labels for genes (variables) in \code{multiExpr}.
}
  \item{MEs}{
If module eigengenes have been calculated before, the user can save some computational time
by inputting them. \code{MEs} should have the same format as \code{multiExpr}.
If they are not given, they will be calculated.
}
  \item{unassdColor}{The label (value in \code{labels})
that represents unassigned genes. Module of this label will
not enter the module eigengene clustering and will not be merged with other modules.}

  \item{impute}{Should missing values be imputed in eigengene calculation? If imputation is disabled, the
presence of \code{NA} entries will cause the eigengene calculation to fail and eigengenes will be
replaced by their hubgene approximation. See \code{\link{moduleEigengenes}} for more details.}

  \item{networkOptions}{
A single list of class \code{\link{NetworkOptions}} giving options for network calculation for all of the
networks, or a \code{\link{multiData}} structure containing one such list for each input data set.
}

  \item{consensusTree}{
A list specifying the consensus calculation. See \code{\link{newConsensusTree}} for details.
}

  \item{calibrateMESimilarities}{
Logical: should module eigengene similarities be calibrated? This setting overrides the calibration options
in \code{consensusTree}.
}
  \item{cutHeight}{
 Maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging. 
}

  \item{iterate}{Controls whether the merging procedure should be repeated until there is no change. If
FALSE, only one iteration will be executed.}

  \item{relabel}{Controls whether, after merging, color labels should be ordered by module size.}

  \item{colorSeq}{Color labels to be used for relabeling. Defaults to the standard color order used
in this package if \code{colors} are not numeric, and to integers starting from 1 if
\code{colors} is numeric.}

  \item{getNewMEs}{Controls whether module eigengenes of merged modules should be calculated and
returned.}
  \item{getNewUnassdME}{When doing module eigengene manipulations, the function does not normally
calculate the eigengene of the 'module' of unassigned ('grey') genes. Setting this option to
\code{TRUE} will force the calculation of the unassigned eigengene in the returned newMEs, but not
in the returned oldMEs.}

  \item{trapErrors}{Controls whether computational errors in calculating module eigengenes, their
dissimilarity, and merging trees should be trapped. If \code{TRUE}, errors will be trapped and the
function will return the input colors. If \code{FALSE}, errors will cause the function to stop.}

  \item{verbose}{Controls verbosity of printed progress messages. 0 means silent, up to (about) 5 the
verbosity gradually increases.}

  \item{indent}{A single non-negative integer controlling indentation of printed messages. 0 means no
indentation, each unit above that adds two spaces. }
}
\details{
  This function merges input modules 
that are closely related. The similarities are quantified by correlations of module eigengenes; a
``consensus'' similarity is calculated using \code{hierarchicalConsensusMEDissimilarity} 
according to the recipe in \code{consensusTree}. Once the
(dis-)similarities are calculated, average linkage hierarchical clustering of the module eigengenes is
performed, the dendrogram is cut at the height \code{cutHeight} and modules on each branch are merged.
The process is (optionally) repeated until no more modules are merged.

If, for a particular module, the module eigengene calculation fails, a hubgene approximation will be
used.

The user should be aware that if a computational error occurs and \code{trapErrors==TRUE},
the returned list (see below) will not contain all of the components returned upon normal execution.
}
\value{
  If no errors occurred, a list with components
  \item{labels}{Labels for the genes corresponding to merged modules. The function attempts to
mimic the mode of the input \code{labels}: if the input \code{labels} is numeric, character and
factor, respectively, so is the output. Note, however, that if the function performs relabeling, a
standard sequence of labels will be used: integers starting at 1 if the input \code{labels} is
numeric, and a sequence of color labels otherwise (see \code{colorSeq} above).}

  \item{dendro}{Hierarchical clustering dendrogram (average linkage) of the eigengenes of the most
recently computed tree. If \code{iterate} was set TRUE, this will be the dendrogram of the merged
modules, otherwise it will be the dendrogram of the original modules.}

  \item{oldDendro}{Hierarchical clustering dendrogram (average linkage) of the eigengenes of the original
modules.}

  \item{cutHeight}{The input cutHeight.}

  \item{oldMEs}{Module eigengenes of the original modules in the sets given by \code{useSets}.}

  \item{newMEs}{Module eigengenes of the merged modules in the sets given by \code{useSets}.}

  \item{allOK}{A logical set to \code{TRUE}.}

  If an error occurred and \code{trapErrors==TRUE}, the list only contains these components:

  \item{colors}{A copy of the input colors.}

  \item{allOK}{a logical set to \code{FALSE}.}
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{multiSetMEs}} for calculation of (consensus) module eigengenes across multiple data sets;

\code{\link{newConsensusTree}} for information about consensus trees;

\code{\link{hierarchicalConsensusMEDissimilarity}} for calculation of hierarchical consensus eigengene
dissimilarity.

}
\keyword{misc}% __ONLY ONE__ keyword per line
