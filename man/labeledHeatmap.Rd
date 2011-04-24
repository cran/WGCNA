\name{labeledHeatmap}
\alias{labeledHeatmap}
\title{ Produce a labeled heatmap plot }
\description{
  Plots a heatmap plot with color legend, row and column annotation, and optional text within th heatmap.
}
\usage{
labeledHeatmap(
  Matrix, 
  xLabels, yLabels = NULL, 
  xSymbols = NULL, ySymbols = NULL, 
  colorLabels = NULL, 
  xColorLabels = FALSE, yColorLabels = FALSE, 
  checkColorsValid = TRUE,
  invertColors = FALSE, 
  setStdMargins = TRUE, 
  xLabelsPosition = "bottom",
  xLabelsAngle = 45,
  xLabelsAdj = 1,
  xColorWidth = 0.05,
  yColorWidth = 0.05,
  colors = NULL, 
  textMatrix = NULL, 
  cex.text = NULL, cex.lab = NULL, 
  cex.lab.x = cex.lab,
  cex.lab.y = cex.lab,
  colors.lab.x = 1,
  colors.lab.y = 1,
  plotLegend = TRUE, ...)
}
\arguments{
  \item{Matrix}{ numerical matrix to be plotted in the heatmap. }

  \item{xLabels}{ labels for the columns. See Details. }

  \item{yLabels}{ labels for the rows. See Details. }

  \item{xSymbols}{ additional labels used when \code{xLabels} are interpreted as colors. See Details. }

  \item{ySymbols}{ additional labels used when \code{yLabels} are interpreted as colors. See Details. }

  \item{colorLabels}{ logical: should \code{xLabels} and \code{yLabels} be interpreted as colors? If
given, overrides \code{xColorLabels} and \code{yColorLabels} below.}

  \item{xColorLabels}{ logical: should \code{xLabels} be interpreted as colors? }

  \item{yColorLabels}{ logical: should \code{yLabels} be interpreted as colors? }

  \item{checkColorsValid}{ logical: should given colors be checked for validity 
against the output of \code{colors()} ? If this argument is \code{FALSE}, invalid color specification
will trigger an error.}

  \item{invertColors}{ logical: should the color order be inverted? }

  \item{setStdMargins}{ logical: should standard margins be set before calling the plot function?
Standard margins depend on \code{colorLabels}: they are wider for text labels and narrower for color
labels. The defaults are static, that is the function does not attempt to guess the optimal margins. }

  \item{xLabelsPosition}{ a character string specifying the position of labels for the columns.
Recognized values are (unique abbreviations of) \code{"top", "bottom"}. }

  \item{xLabelsAngle}{ angle by which the column labels should be rotated. }

  \item{xLabelsAdj}{ justification parameter for column labels. See \code{\link{par}} and the
description of parameter \code{"adj"}. }

  \item{xColorWidth}{ width of the color labels for the x axis expressed as a fraction of the smaller of
the range of the x and y axis. }

  \item{yColorWidth}{ width of the color labels for the y axis expressed as a fraction of the smaller of
the range of the x and y axis. }

  \item{colors}{ color pallette to be used in the heatmap. Defaults to \code{\link{heat.colors}}. }

  \item{textMatrix}{ optional matrix of text entries of the same dimensions as \code{Matrix}. }

  \item{cex.text}{ character expansion factor for \code{textMatrix}. }

  \item{cex.lab}{ character expansion factor for text labels labeling the axes. }

  \item{cex.lab.x}{ character expansion factor for text labels labeling the x axis. Overrides \code{cex.lab}
above. }

  \item{cex.lab.y}{ character expansion factor for text labels labeling the y axis. Overrides \code{cex.lab}  
above. }

  \item{colors.lab.x}{colors for character labels or symbols along x axis.} 

  \item{colors.lab.y}{colors for character labels or symbols along y axis.} 

  \item{plotLegend}{ logical: should a color legend be plotted? }

  \item{\dots}{ other arguments to function \code{\link{heatmap}}. }
}
\details{
  The function basically plots a standard heatmap plot of the given \code{Matrix} and embellishes it with
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
}
\value{
  None.
}
\author{ Peter Langfelder}
\seealso{ \code{\link{heatmap}}, \code{\link{colors}} }

\examples{

# This example illustrates 4 main ways of annotating columns and rows of a heatmap.
# Copy and paste the whole example into an R session with an interactive plot window;
# alternatively, you may replace the command sizeGrWindow below by opening another graphical device such
# as pdf.

# Generate a matrix to be plotted

nCol = 8; nRow = 7;
mat = matrix(runif(nCol*nRow, min = -1, max = 1), nRow, nCol);

rowColors = standardColors(nRow);
colColors = standardColors(nRow + nCol)[(nRow+1):(nRow + nCol)];

rowColors;
colColors;

sizeGrWindow(9,7)
par(mfrow = c(2,2))
par(mar = c(4, 5, 4, 6));

# Label rows and columns by text:

labeledHeatmap(mat, xLabels = colColors, yLabels = rowColors, 
               colors = greenWhiteRed(50),
               setStdMargins = FALSE, 
               textMatrix = signif(mat, 2),
               main = "Text-labeled heatmap");

# Label rows and columns by colors:

rowLabels = paste("ME", rowColors, sep="");
colLabels = paste("ME", colColors, sep="");

labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
               colorLabels = TRUE,
               colors = greenWhiteRed(50),
               setStdMargins = FALSE,
               textMatrix = signif(mat, 2),
               main = "Color-labeled heatmap");

# Mix text and color labels:

rowLabels[3] = "Row 3";
colLabels[1] = "Column 1";

labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
               colorLabels = TRUE,
               colors = greenWhiteRed(50),
               setStdMargins = FALSE,
               textMatrix = signif(mat, 2), 
               main = "Mix-labeled heatmap");

# Color labels and additional text labels

rowLabels = paste("ME", rowColors, sep="");
colLabels = paste("ME", colColors, sep="");

extraRowLabels = paste("Row", c(1:nRow));
extraColLabels = paste("Column", c(1:nCol));

# Extend margins to fit all labels
par(mar = c(6, 6, 4, 6));
labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
               xSymbols = extraColLabels,
               ySymbols = extraRowLabels,
               colorLabels = TRUE,
               colors = greenWhiteRed(50),
               setStdMargins = FALSE,
               textMatrix = signif(mat, 2),
               main = "Text- + color-labeled heatmap");

}

\keyword{ hplot }% __ONLY ONE__ keyword per line
