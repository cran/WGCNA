\name{plotColorUnderTree}
\alias{plotColorUnderTree}
\alias{plotOrderedColors}
\title{Plot color rows in a given order, for example under a dendrogram}
\description{
  Plot color rows encoding information about objects in a given order, for example the order of 
a clustering dendrogram, 
usually below the dendrogram or a barplot.
}
\usage{
plotOrderedColors(
   order, 
   colors, 
   rowLabels = NULL, 
   rowWidths = NULL, 
   rowText = NULL,
   rowTextAlignment = c("left", "center", "right"),
   rowTextIgnore = NULL,
   textPositions = NULL, 
   addTextGuide = TRUE,
   cex.rowLabels = 1, 
   cex.rowText = 0.8,
   startAt = 0,
   ...)

plotColorUnderTree(
   dendro, 
   colors,
   rowLabels = NULL,
   rowWidths = NULL,
   rowText = NULL,
   rowTextAlignment = c("left", "center", "right"),
   rowTextIgnore = NULL,
   textPositions = NULL,
   addTextGuide = TRUE,
   cex.rowLabels = 1,
   cex.rowText = 0.8,
   ...)

}
\arguments{

  \item{order}{A vector giving the order of the objects. Must have the same length as \code{colors} if
\code{colors} is a vector, or as the number of rows if \code{colors} is a matrix or data frame.}

  \item{dendro}{A hierarchical clustering dendrogram such one returned by \code{\link{hclust}}.}

  \item{colors}{Coloring of objects on the dendrogram. Either a vector (one color per object) or a
matrix (can also be an array or a data frame)
with each column giving one color per object. Each column will be plotted as a horizontal row of colors
under the dendrogram.}

  \item{rowLabels}{Labels for the colorings given in \code{colors}. The labels will be printed to the
left of the color rows in the plot. If the argument is given, it must be a vector of length
equal to the number of columns in \code{colors}. If not given, \code{names(colors)}
will be used if available. If not, sequential numbers
starting from 1 will be used.}

  \item{rowWidths}{ Optional specification of relative row widths for the color and text (if given) rows.
Need not sum to 1. }

  \item{rowText}{Optional labels to identify colors in the color rows. 
If given, must be of the same dimensions as \code{colors}. Each label that occurs will be displayed
once.}

  \item{rowTextAlignment}{Character string specifying whether the labels should be left-justified to the
start of the largest block of each label, centered in the middle, or right-justified to the end of the
largest block.}

  \item{rowTextIgnore}{Optional specifications of labels that should be ignored when displaying them using
\code{rowText} above. }

  \item{textPositions}{optional numeric vector of the same length as the number of columns in \code{rowText}
giving the color rows under which the text rows should appear.}

  \item{addTextGuide}{ logical: should guide lines be added for the text rows (if given)? }

  \item{cex.rowLabels}{Font size scale factor for the row labels. See \code{\link[graphics]{par}}.}

  \item{cex.rowText}{ character expansion factor for text rows (if given). }

  \item{startAt}{A numeric value indicating where in relationship to the left edge of the plot the center
of the first rectangle should be. Useful values are 0 if ploting color under a dendrogram, and 0.5 if
ploting colors under a barplot. }

  \item{\dots}{Other parameters to be passed on to the plotting method (such as \code{main} for the main
title etc).}
}
\details{
    It is often useful to plot dendrograms or other plots (e.g., barplots) 
    of objects together with additional information about the
objects, for example module assignment (by color) that was obtained by cutting a hierarchical
dendrogram or external color-coded measures such as gene significance.
This function provides a way to do so. The calling code should section the screen into two 
(or more) parts, plot the dendrogram (via \code{plot(hclust)}) or other information 
in the upper section and use this function
to plot color annotation in the order corresponding to the dendrogram in the lower section. 
}
\value{ None.
}
\note{ This function replaces \code{plotHclustColors} in package \code{moduleColor}.
}
\author{
  Steve Horvath \email{SHorvath@mednet.ucla.edu} and Peter Langfelder \email{Peter.Langfelder@gmail.com}
}
\seealso{\code{\link[dynamicTreeCut]{cutreeDynamic}} for module detection in a dendrogram;

\code{\link{plotDendroAndColors}} for automated plotting of dendrograms and colors in one step.}

\keyword{hplot}

