\name{automaticNetworkScreening}
\alias{automaticNetworkScreening}
\title{ One-step automatic network gene screening }
\description{
   This function performs gene screening based on a given trait and gene network properties
}
\usage{
automaticNetworkScreening(
   datExpr, 
   y, 
   power = 6, 
   networkType = "unsigned", 
   detectCutHeight = 0.995,
   minModuleSize = min(20, ncol(as.matrix(datExpr))/2), 
   datME = NULL, 
   getQValues = TRUE, 
   ...)
}
\arguments{
  \item{datExpr}{ data frame containing the expression data, columns corresponding to genes and rows to
samples }

  \item{y}{ vector containing trait values for all samples in \code{datExpr}}

  \item{power}{ soft thresholding power used in network construction }

  \item{networkType}{ character string specifying network type. Allowed values are (unique abbreviations
of) \code{"unsigned"}, \code{"signed"}, \code{"hybrid"}. }

  \item{detectCutHeight}{ cut height of the gene hierarchical clustering dendrogram. See
\code{cutreeDynamic} for details. }
  \item{minModuleSize}{ minimum module size to be used in module detection procedure. }
  \item{datME}{ optional specification of module eigengenes. A data frame whose columns are the module
eigengenes. If given, module analysis will not be performed. }
  \item{getQValues}{logical: should q-values (local FDR) be calculated? }
  \item{...}{other arguments to the module identification function \code{\link{blockwiseModules}}}
}
\details{
  Network screening is a method for identifying genes that have a high gene significance and are members
of important modules at the same time.
  If \code{datME} is given, the function calls \code{\link{networkScreening}} with the default
parameters. If \code{datME} is not given, module eigengenes are first calculated using network analysis
based on supplied parameters.
}

\value{
  A list with the following components:
  \item{networkScreening}{a data frame containing results of the network screening procedure. See
\code{\link{networkScreening}} for more details.}
  \item{datME}{ calculated module eigengenes (or a copy of the input \code{datME}, if given).}
  \item{hubGeneSignificance}{ hub gene significance for all calculated modules. See
\code{\link{hubGeneSignificance}}. }
}

\author{ Steve Horvath }

\seealso{ 
\code{\link{networkScreening}}, \code{\link{hubGeneSignificance}},
\code{\link{networkScreening}}, \code{\link[dynamicTreeCut]{cutreeDynamic}} }

\keyword{ misc }% __ONLY ONE__ keyword per line

