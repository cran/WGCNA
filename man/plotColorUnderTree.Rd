\name{plotColorUnderTree}
\alias{plotColorUnderTree}
\title{Plot color rows under a dendrogram}
\description{
  Plot color rows encoding information about objects in a clustering dendrogram, 
usually below the dendrogram.
}
\usage{
plotColorUnderTree(
   dendro, 
   colors, 
   rowLabels = NULL, 
   rowWidths = NULL, 
   rowText = NULL,
   textPositions = NULL, 
   addTextGuide = TRUE,
   cex.rowLabels = 1, 
   cex.rowText = 0.8,
   ...)
}
\arguments{

  \item{dendro}{A dendrogram such as returned by \code{\link[stats]{hclust}}. }

  \item{colors}{Coloring of objects on the dendrogram. Either a vector (one color per object) or a
matrix (can also be an array or a data frame)
with each column giving one color per object. Each column will be plotted as a horizontal row of colors
under the dendrogram.}

  \item{rowLabels}{Labels for the colorings given in \code{colors}. The labels will be printed to the
left of the color rows in the plot. If the argument is given, it must be a vector of length
equal to the number of columns in \code{colors}. If not given, \code{names(colors)}
will be used if available. If not, sequential numbers
starting from 1 will be used.}

  \item{rowWidths}{ optional specification of relative row widths for the color and text (if given) rows.
Need not sum to 1. }

  \item{rowText}{Optional labels to identify colors in the color rows. 
If given, must be of the same dimensions as \code{colors}. Each label that occurs will be displayed
once.}

  \item{textPositions}{optional numeric vector of the same length as the number of columns in \code{rowText}
giving the color rows under which the text rows should appear.}

  \item{addTextGuide}{ logical: should guide lines be added for the text rows (if given)? }

  \item{cex.rowLabels}{Font size scale factor for the row labels. See \code{\link[graphics]{par}}.}

  \item{cex.rowText}{ character expansion factor for text rows (if given). }

  \item{\dots}{Other parameters to be passed on to the plotting method (such as \code{main} for the main
title etc).}
}
\details{
    It is often useful to plot dendrograms of objects together with additional information about the
objects, for example module assignment (by color) that was obtained by cutting a hierarchical
dendrogram or external color-coded measures such as gene significance.
This function provides a way to do so. The calling code should section the screen into two 
(or more) parts, plot the dendrogram (via \code{plot(hclust)}) in the upper section and use this function
to plot color annotation in the order corresponding to the dendrogram in the lower section. 
}
\value{ None.
}
\note{ This function is identical to the function \code{plotHclustColors} in package \code{moduleColor}.
}
\author{
  Steve Horvath \email{SHorvath@mednet.ucla.edu} and Peter Langfelder \email{Peter.Langfelder@gmail.com}
}
\seealso{\code{\link[dynamicTreeCut]{cutreeDynamic}} for module detection in a dendrogram;

\code{\link{plotDendroAndColors}} for automated plotting of dendrograms and colors in one step.}

\keyword{hplot}

