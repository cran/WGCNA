\name{mergeCloseModules}
\alias{mergeCloseModules}
\title{Merge close modules in gene expression data}
\description{
  Merges modules in gene expression networks that are too close as measured by the correlation of their
eigengenes.
}
\usage{
mergeCloseModules(
 # input data
  exprData, colors,

  # Optional starting eigengenes
  MEs = NULL,

  # Optional restriction to a subset of all sets
  useSets = NULL,

  # If missing data are present, impute them?
  impute = TRUE,

  # Input handling options
  checkDataFormat = TRUE,
  unassdColor = if (is.numeric(colors)) 0 else "grey",

  # Options for eigengene network construction
  corFnc = cor, corOptions = list(use = 'p'),
  useAbs = FALSE,

  # Options for constructing the consensus
  equalizeQuantiles = FALSE,
  quantileSummary = "mean",
  consensusQuantile = 0,

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
\arguments{
  \item{exprData}{Expression data, either a single data frame with rows corresponding to samples and
columns to genes, or in a multi-set format (see \code{\link{checkSets}}). See
\code{checkDataStructure} below.  }

  \item{colors}{A vector (numeric, character or a factor) giving module colors for genes. 
The method only makes sense when genes have the
same color label in all sets, hence a single vector. }

  \item{MEs}{If module eigengenes have been calculated before, the user can save some computational time
by inputting them. \code{MEs} should have the same format as \code{exprData}. 
If they are not given, they will be calculated.}

  \item{useSets}{A vector of scalar allowing the user to specify which sets will be used to calculate the
consensus dissimilarity of module eigengenes. Defaults to all given sets. } 

  \item{impute}{Should missing values be imputed in eigengene calculation? If imputation is disabled, the
presence of \code{NA} entries will cause the eigengene calculation to fail and eigengenes will be
replaced by their hubgene approximation. See \code{\link{moduleEigengenes}} for more details.}

  \item{checkDataFormat}{If TRUE, the function will check \code{exprData} and \code{MEs} for correct
multi-set structure. If single set data is given, it will be converted into a format usable for the
function. If FALSE, incorrect structure of input data will trigger an error.}

  \item{unassdColor}{Specifies the string that labels unassigned genes. Module of this color will
not enter the module eigengene clustering and will not be merged with other modules.}

  \item{corFnc}{Correlation function to be used to calculate correlation of module eigengenes. }

  \item{corOptions}{Can be used to specify options to the correlation function, in addition to argument
\code{x} which is used to pass the actual data to calculate the correlation of.}

  \item{useAbs}{Specifies whether absolute value of correlation or plain correlation (of module
eigengenes) should be used in calculating module dissimilarity.}

  \item{equalizeQuantiles}{Logical: should quantiles of the eigengene dissimilarity matrix be equalized
("quantile normalized")? The default is \code{FALSE} for reproducibility of old code; when there are many
eigengenes (e.g., at least 50), better results may be achieved if quantile equalization is used.}

  \item{quantileSummary}{One of \code{"mean"} or \code{"median"}. Controls how a reference dissimilarity
is computed from the input ones (using mean or median, respectively).} 

  \item{consensusQuantile}{A number giving the desired quantile to use in the consensus similarity
calculation (see details).}

  \item{cutHeight}{Maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging.}

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
that are closely related. The similarities are measured by correlations of module eigengenes; a
``consensus'' measure is defined as the ``consensus quantile'' 
over the corresponding relationship in each set. Once the
(dis-)similarity is calculated, average linkage hierarchical clustering of the module eigengenes is
performed, the dendrogram is cut at the height \code{cutHeight} and modules on each branch are merged.
The process is (optionally) repeated until no more modules are merged.

If, for a particular module, the module eigengene calculation fails, a hubgene approximation will be
used. 

The user should be aware that if a computational error occurs and \code{trapErrors==TRUE}, 
the returned list (see below) will not contain all of the components returned upon normal execution.
}
\value{
  If no errors occurred, a list with components
  \item{colors}{Color labels for the genes corresponding to merged modules. The function attempts to
mimic the mode of the input \code{colors}: if the input \code{colors} is numeric, character and
factor, respectively, so is the output. Note, however, that if the fnction performs relabeling, a
standard sequence of labels will be used: integers starting at 1 if the input \code{colors} is
numeric, and a sequence of color labels otherwise (see \code{colorSeq} above).}

  \item{dendro}{Hierarchical clustering dendrogram (average linkage) of the eigengenes of the most
recently computed tree. If \code{iterate} was set TRUE, this will be the dendrogram of the merged
modules, otherwise it will be the dendrogram of the original modules.}

  \item{oldDendro}{Hierarchical clustering dendrogram (average linkage) of the eigengenes of the original
modules.}

  \item{cutHeight}{The input cutHeight.}

  \item{oldMEs}{Module eigengenes of the original modules in the sets given by \code{useSets}.}

  \item{newMEs}{Module eigengenes of the merged modules in the sets given by \code{useSets}.}

  \item{allOK}{A boolean set to \code{TRUE}.}

  If an error occurred and \code{trapErrors==TRUE}, the list only contains these components:

  \item{colors}{A copy of the input colors.}

  \item{allOK}{a boolean set to \code{FALSE}.}
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
%\seealso{\code{\link{moduleEigengenes}}, \code{\link{multiSetMEs}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
