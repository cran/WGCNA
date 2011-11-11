\name{multiData.eigengeneSignificance}
\alias{multiData.eigengeneSignificance}

\title{
Eigengene significance across multiple sets
}
\description{
This function calculates eigengene significance and the associated significance statistics (p-values,
q-values etc) across several data sets.
}
\usage{
multiData.eigengeneSignificance(
  multiData, multiTrait, 
  moduleLabels, multiEigengenes = NULL, 
  useModules = NULL, 
  corAndPvalueFnc = corAndPvalue, corOptions = list(), 
  corComponent = "cor", 
  getQvalues = FALSE, setNames = NULL, 
  excludeGrey = TRUE, greyLabel = ifelse(is.numeric(moduleLabels), 0, "grey"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiData}{
   Expression data (or other data) in multi-set format (see \code{\link{checkSets}}). A vector of lists; in
each list there must be a component named \code{data} whose content
is a matrix or dataframe or array of dimension 2.
}
  \item{multiTrait}{
   Trait or ourcome data in multi-set format. Only one trait is allowed; consequesntly, the \code{data}
component of each component list can be either a vector or a data frame (matrix, array of dimension 2).
}
  \item{moduleLabels}{
  Module labels: one label for each gene in \code{multiExpr}.
}
  \item{multiEigengenes}{
  Optional eigengenes of modules specified in \code{moduleLabels}. If not given, will be calculated from
\code{multiExpr}.
}
  \item{useModules}{
  Optional specification of module labels to which the analysis should be restricted. This could be useful
if there are many modules, most of which are not interesting. Note that the "grey" module cannot be used
with \code{useModules}.}

  \item{corAndPvalueFnc}{
  Function that calculates associations between expression profiles and eigengenes. See details.
}
  \item{corOptions}{
  List giving additional arguments to function \code{corAndPvalueFnc}. See details.
}
  \item{corComponent}{
  Name of the component of output of \code{corAndPvalueFnc} that contains the actual correlation.
}
  \item{getQvalues}{
   logical: should q-values (estimates of FDR) be calculated?
}
  \item{setNames}{
   names for the input sets. If not given, will be taken from \code{names(multiExpr)}. If those are
\code{NULL} as well, the names will be \code{"Set_1", "Set_2", ...}.
}
   \item{excludeGrey}{
   logical: should the grey module be excluded from the kME tables? Since the grey module is typically not a
real module, it makes little sense to report kME values for it.
}
   \item{greyLabel}{
   label that labels the grey module.
}
}
\details{
This is a convenience function that calculates module eigengene significances (i.e., correlations of module
eigengenes with a given trait) across all sets in a multi-set analysis. Also returned are p-values, Z
scores, numbers of present (i.e., non-missing) observations for each significance, and optionally the
q-values (false discovery rates) corresponding to the p-values.

The function \code{corAndPvalueFnc} is currently
is expected to accept arguments \code{x} (gene expression profiles) and \code{y} (eigengene expression
profiles).  Any additional arguments can be passed via \code{corOptions}. 

The function \code{corAndPvalueFnc} should return a list which at the least contains (1) a matrix
of associations of genes and eigengenes (this component should have the name given by \code{corComponent}), 
and (2) a matrix of the corresponding p-values, named "p" or "p.value". Other components are optional but
for full functionality should include
(3) \code{nObs} giving the number of observations for each association (which is the number of samples less
number of missing data - this can in principle vary from association to association), and (4) \code{Z}
giving a Z static for each observation. If these are missing, \code{nObs} is calculated in the main
function, and calculations using the Z statistic are skipped.
}
\value{
A list containing the following components. Each component is a matrix in which the rows correspond to
module eigengenes and columns to data sets. Row and column names are set appropriately.
\item{eigengeneSignificance}{Module eigengene significance.}
\item{p.value}{p-values (returned by \code{corAndPvalueFnc}). }
\item{q.value}{q-values corresponding to the p-values above. Only returned in input \code{getWvalues} is
\code{TRUE}. }
\item{Z}{Z statistics (if returned by \code{corAndPvalueFnc}). }
\item{nObservations}{Number of non-missing observations in each correlation/p-value.}
}
\author{
Peter Langfelder
}
\keyword{misc}
