\name{labeledHeatmap.multiPage}
\alias{labeledHeatmap.multiPage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Labeled heatmap divided into several separate plots.
}
\description{
This function produces labaled heatmaps divided into several plots. This is useful for large heatmaps where
labels on individual columns and rows may become unreadably small (or overlap).
}
\usage{
labeledHeatmap.multiPage(
   # Input data and ornaments
   Matrix, 
   xLabels, yLabels = NULL,
   xSymbols = NULL, ySymbols = NULL,
   textMatrix = NULL, 

   # Paging options
   rowsPerPage = NULL, maxRowsPerPage = 20, 
   colsPerPage = NULL, maxColsPerPage = 10, 
   addPageNumberToMain = TRUE, 

   # Further arguments to labeledHeatmap
   zlim = NULL,
   signed = TRUE, 
   main = "", 

   # Separator line specification                   
   verticalSeparator.x = NULL,
   verticalSeparator.col = 1,
   verticalSeparator.lty = 1,
   verticalSeparator.lwd = 1,
   verticalSeparator.ext = 0,

   horizontalSeparator.y = NULL,
   horizontalSeparator.col = 1,
   horizontalSeparator.lty = 1,
   horizontalSeparator.lwd = 1,
   horizontalSeparator.ext = 0,

   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Matrix}{ numerical matrix to be plotted in the heatmap. }

  \item{xLabels}{ labels for the columns. See Details. }

  \item{yLabels}{ labels for the rows. See Details. }

  \item{xSymbols}{ additional labels used when \code{xLabels} are interpreted as colors. See Details. }

  \item{ySymbols}{ additional labels used when \code{yLabels} are interpreted as colors. See Details. }

  \item{textMatrix}{ optional text entries for each cell. Either a matrix of the same dimensions as 
\code{Matrix} or a vector of the same length as the number of entries in \code{Matrix}. }


  \item{rowsPerPage}{ optional list in which each component is a vector specifying which rows should appear
together in each plot. If not given, will be generated automatically based on \code{maxRowsPerPage} below
and the number of rows in \code{Matrix}.
}
  \item{maxRowsPerPage}{ integer giving maximum number of rows appearing on each plot (page).  
}
  \item{colsPerPage}{ optional list in which each component is a vector specifying which columns should appear 
together in each plot. If not given, will be generated automatically based on \code{maxColsPerPage} below
and the number of rows in \code{Matrix}.

}
  \item{maxColsPerPage}{
integer giving maximum number of columns appearing on each plot (page).
}
  \item{addPageNumberToMain}{
logical: should plot/page number be added to the \code{main} title of each plot?
}
  \item{zlim}{
Optional specification of the extreme values for the color scale. If not given, will be determined from the
input \code{Matrix}.
}
  \item{main}{
Main title for each plot/page, optionally with the plot/page number added.
}
  \item{signed}{
logical: should the input \code{Matrix} be converted to colors using a scale centered at zero?
}
  \item{verticalSeparator.x}{indices of columns after which separator lines (vertical lines between columns)
      should be drawn. \code{NULL} means no lines will be drawn.}
  \item{verticalSeparator.col}{color(s) of the vertical separator lines. Recycled if need be. }
  \item{verticalSeparator.lty}{line type of the vertical separator lines. Recycled if need be. }
  \item{verticalSeparator.lwd}{line width of the vertical separator lines. Recycled if need be. }
  \item{verticalSeparator.ext}{number giving the extension of the separator line into the margin as a fraction
of the margin width. 0 means no extension, 1 means extend all the way through the margin. }

  \item{horizontalSeparator.y}{indices of columns after which separator lines (horizontal lines between
columns)
      should be drawn. \code{NULL} means no lines will be drawn.}
  \item{horizontalSeparator.col}{ color(s) of the horizontal separator lines. Recycled if need be. }
  \item{horizontalSeparator.lty}{line type of the horizontal separator lines. Recycled if need be. }
  \item{horizontalSeparator.lwd}{line width of the horizontal separator lines. Recycled if need be. }
  \item{horizontalSeparator.ext}{number giving the extension of the separator line into the margin as a
fraction
of the margin width. 0 means no extension, 1 means extend all the way through the margin. }

  \item{\dots}{
other arguments to function \code{\link{labeledHeatmap}}.
}
}
\details{
The function \code{\link{labeledHeatmap}} is used to produce each plot/page; most arguments are described
in more detail in the help file for that function.

In each plot/page \code{\link{labeledHeatmap}} plots a standard heatmap plot of an appropriate 
sub-rectangle of \code{Matrix} and embellishes it with
row and column labels and/or with text within the heatmap entries. Row and column labels can be either
character strings or color squares, or both.

To get simple text labels, use \code{colorLabels=FALSE} and pass the desired row and column labels in
\code{yLabels} and \code{xLabels}, respectively.

To label rows and columns by color squares, use
\code{colorLabels=TRUE}; \code{yLabels} and \code{xLabels} are then expected to represent valid colors.
For reasons of compatibility with other functions, each entry in \code{yLabels} and \code{xLabels} is
expected to consist of a color designation preceded by 2 characters: an example would be
\code{MEturquoise}. The first two characters can be arbitrary, they are stripped.
Any labels that do not represent valid colors will be considered text labels and printed in full,
allowing the user to mix text and color labels.

It is also possible to label rows and columns by both color squares and additional text annotation. To
achieve this, use the above technique to get color labels and, additionally, pass the desired text
annotation in the \code{xSymbols} and \code{ySymbols} arguments.


If \code{rowsPerPage} (\code{colsPerPage}) is not given, rows (columns) are allocated automatically as
uniformly as possible, in contiguous blocks of size at most \code{maxRowsPerPage} (\code{maxColsPerPage}). 
The allocation is performed by the function \code{\link{allocateJobs}}.
}
\value{
  None.
}
\author{
Peter Langfelder
}

\seealso{
The workhorse function \code{\link{labeledHeatmap}} for the actual heatmap plot;

function \code{\link{allocateJobs}} for the allocation of rows/columns to each plot.

}
% R documentation directory.
\keyword{ misc }
