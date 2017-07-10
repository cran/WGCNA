\name{plotDendroAndColors}
\alias{plotDendroAndColors}
\title{ Dendrogram plot with color annotation of objects }
\description{
  This function plots a hierarchical clustering dendrogram and color annotation(s) of objects in the
dendrogram underneath.
}
\usage{
plotDendroAndColors(
  dendro, 
  colors, 
  groupLabels = NULL, 
  rowText = NULL,
  rowTextAlignment = c("left", "center", "right"),
  rowTextIgnore = NULL,
  textPositions = NULL,
  setLayout = TRUE, 
  autoColorHeight = TRUE, 
  colorHeight = 0.2, 
  colorHeightBase = 0.2, 
  colorHeightMax = 0.6,
  rowWidths = NULL, 
  dendroLabels = NULL, 
  addGuide = FALSE, guideAll = FALSE, 
  guideCount = 50, guideHang = 0.2, 
  addTextGuide = FALSE,
  cex.colorLabels = 0.8, cex.dendroLabels = 0.9, 
  cex.rowText = 0.8,
  marAll = c(1, 5, 3, 1), saveMar = TRUE, 
  abHeight = NULL, abCol = "red", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendro}{ a hierarchical clustering dendrogram such as one produced by
\code{\link[stats]{hclust}}. }

  \item{colors}{ Coloring of objects on the dendrogram. Either a vector (one color per object) or a
matrix (can also be an array or a data frame)
with each column giving one color per object. Each column will be plotted as a horizontal row of colors
under the dendrogram. }

  \item{groupLabels}{ Labels for the colorings given in \code{colors}. The labels will be printed to the
left of the color rows in the plot. If the argument is given, it must be a vector of length
equal to the number of columns in \code{colors}. If not given, \code{names(colors)}
will be used if available. If not, sequential numbers starting from 1 will be used.}

  \item{rowText}{Optional labels to identify colors in the color rows. 
If given, must be either the same dimensions as \code{colors} or must have the same number of rows and
\code{textPositions} must be used to specify which columns of \code{colors} each column of \code{rowText}
corresponds to. Each label that occurs will be displayed 
once, under the largest continuous block of the corresponding \code{colors}.}

  \item{rowTextAlignment}{Character string specifying whether the labels should be left-justified to the
start of the largest block of each label, centered in the middle, or right-justified to the end of the 
largest block.}

  \item{rowTextIgnore}{Optional specifications of labels that should be ignored when displaying them using 
\code{rowText} above. }

  \item{textPositions}{optional numeric vector of the same length as the number of columns in \code{rowText}
giving the color rows under which the text rows should appear.} 

  \item{setLayout}{ logical: should the plotting device be partitioned into a standard layout?
If \code{FALSE}, the user is responsible for partitioning. The function expects two regions of the same
width, the first one immediately above the second one. }

  \item{autoColorHeight}{ logical: should the height of the color area below the dendrogram be
automatically adjusted for the number of traits? Only effective if \code{setLayout} is \code{TRUE}. }

  \item{colorHeight}{ specifies the height of the color area under dendrogram as a fraction of the
height of the dendrogram area. Only effective when \code{autoColorHeight} above is \code{FALSE}. }

  \item{colorHeightBase}{when \code{autoColorHeight} is \code{TRUE}, this specifies the minimum height of the color
area (the height when there is one color row).}

  \item{colorHeightMax}{when \code{autoColorHeight} is \code{TRUE}, this specifies the maximum height of the color
area (the height when there are many color rows).}

  \item{rowWidths}{ optional specification of relative row widths for the color and text (if given) rows.
Need not sum to 1. }

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

  \item{addTextGuide}{ logical: should guide lines be added for the text rows (if given)? }

  \item{cex.colorLabels}{ character expansion factor for trait labels. }

  \item{cex.dendroLabels}{ character expansion factor for dendrogram (sample) labels. }

  \item{cex.rowText}{ character expansion factor for text rows (if given). }

  \item{marAll}{ a vector of length 4 giving the bottom, left, top and right margins of the combined
plot. There is no margin between the dendrogram and the color plot underneath. }

  \item{saveMar}{ logical: save margins setting before starting the plot and restore on exit? }

  \item{abHeight}{ optional specification of the height for a horizontal line in the dendrogram, see
\code{\link{abline}}. }

  \item{abCol}{ color for plotting the horizontal line. }

  \item{\dots}{ other graphical parameters to \code{\link{plot.hclust}}. }
}
\details{
  The function slits the plotting device into two regions, plots the given dendrogram in the upper
region, then plots color rows in the region below the dendrogram. 
}
\value{
   None.
}
\author{ Peter Langfelder }
\seealso{ \code{\link{plotColorUnderTree}} }
\keyword{ hplot }
