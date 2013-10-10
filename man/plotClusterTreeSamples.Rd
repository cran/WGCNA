\name{plotClusterTreeSamples}
\alias{plotClusterTreeSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Annotated clustering dendrogram of microarray samples }
\description{
  This function plots an annotated clustering dendorgram of microarray samples.
}
\usage{
plotClusterTreeSamples(
  datExpr, 
  y = NULL, 
  traitLabels = NULL, 
  yLabels = NULL,
  main = if (is.null(y)) "Sample dendrogram" else 
                         "Sample dendrogram and trait indicator", 
  setLayout = TRUE, autoColorHeight = TRUE, colorHeight = 0.3,
  dendroLabels = NULL, 
  addGuide = FALSE, guideAll = TRUE, 
  guideCount = NULL, guideHang = 0.2, 
  cex.traitLabels = 0.8, 
  cex.dendroLabels = 0.9, 
  marAll = c(1, 5, 3, 1), 
  saveMar = TRUE, 
  abHeight = NULL, abCol = "red", 
  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{datExpr}{  a data frame containing expression data, with rows corresponding to samples and
columns to genes. Missing values are allowed and will be ignored. }

  \item{y}{ microarray sample trait. Either a vector with one entry per sample, or a matrix in which each
column corresponds to a (different) trait and each row to a sample. }

  \item{traitLabels}{ labels to be printed next to the color rows depicting sample traits. Defaults to
column names of \code{y}.}

  \item{yLabels}{Optional labels to identify colors in the row identifying the sample classes.
If given, must be of the same dimensions as \code{y}. Each label that occurs will be displayed
once.}

  \item{main}{ title for the plot. }

  \item{setLayout}{ logical: should the plotting device be partitioned into a standard layout? 
If \code{FALSE}, the user is responsible for partitioning. The function expects two regions of the same
width, the first one immediately above the second one. }

  \item{autoColorHeight}{ logical: should the height of the color area below the dendrogram be
automatically adjusted for the number of traits? Only effective if \code{setLayout} is \code{TRUE}. }

  \item{colorHeight}{ Specifies the height of the color area under dendrogram as a fraction of the
height of the dendrogram area. Only effective when \code{autoColorHeight} above is \code{FALSE}. }

  \item{dendroLabels}{ dendrogram labels. Set to \code{FALSE} to disable dendrogram labels altogether;
set to \code{NULL} to use row labels of \code{datExpr}. }

  \item{addGuide}{ logical: should vertical "guide lines" be added to the dendrogram plot? The lines make
it easier to identify color codes with individual samples. }

  \item{guideAll}{ logical: add a guide line for every sample? Only effective for \code{addGuide} set
\code{TRUE}. }

  \item{guideCount}{ number of guide lines to be plotted. Only effective when \code{addGuide} is
\code{TRUE} and \code{guideAll} is \code{FALSE}.  }

  \item{guideHang}{ fraction of the dendrogram height to leave between the top end of the guide line and
the dendrogram merge height. If the guide lines overlap with dendrogram labels, increase \code{guideHang}
to leave more space for the labels. }

  \item{cex.traitLabels}{ character expansion factor for trait labels. }

  \item{cex.dendroLabels}{ character expansion factor for dendrogram (sample) labels. }

  \item{marAll}{ a 4-element vector giving the bottom, left, top and right margins around the combined
plot. Note that this is not the same as setting the margins via a call to \code{\link{par}}, because the
bottom margin of the dendrogram and the top margin of the color underneath are always zero. }

  \item{saveMar}{ logical: save margins setting before starting the plot and restore on exit? }

  \item{abHeight}{ optional specification of the height for a horizontal line in the dendrogram, see
\code{\link{abline}}. }

  \item{abCol}{ color for plotting the horizontal line. }

  \item{\dots}{ other graphical parameters to \code{\link{plot.hclust}}. }
}
\details{
  The function generates an average linkage hierarchical clustering dendrogram (see
\code{\link[stats]{hclust}}) of samples from the given expression data, using Eclidean distance of
samples. The dendrogram is plotted together with color annotation for the samples. 

The trait \code{y} must be numeric. If \code{y} is integer, the colors will correspond to values. If
\code{y} is continouos, it will be dichotomized to two classes, below and above median.

}
\value{
  None.
}
\author{ Steve Horvath and Peter Langfelder}
\seealso{ \code{\link[stats]{dist}}, \code{\link[stats]{hclust}}, \code{\link{plotDendroAndColors}} }
\keyword{ hplot }
\keyword{ misc }
