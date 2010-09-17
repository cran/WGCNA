\name{plotEigengeneNetworks}
\alias{plotEigengeneNetworks}
\title{ Eigengene network plot }
\description{
  This function plots dendrogram and eigengene representations of (consensus) eigengenes networks. 
  In the case of conensus eigengene networks the function also plots pairwise preservation  measures
between consensus networks in different sets. 
}
\usage{
plotEigengeneNetworks(
  multiME, 
  setLabels, 
  letterSubPlots = FALSE, Letters = NULL, 
  excludeGrey = TRUE, greyLabel = "grey", 
  plotDendrograms = TRUE, plotHeatmaps = TRUE, 
  setMargins = TRUE, marDendro = NULL, marHeatmap = NULL,
  colorLabels = TRUE, signed = TRUE, 
  heatmapColors = NULL, 
  plotAdjacency = TRUE,
  printAdjacency = FALSE, cex.adjacency = 0.9,
  coloredBarplot = TRUE, barplotMeans = TRUE, barplotErrors = FALSE, 
  plotPreservation = "standard",
  zlimPreservation = c(0, 1), 
  printPreservation = FALSE, cex.preservation = 0.9, 
  ...)
}
\arguments{
  \item{multiME}{ either a single data frame containing the module eigengenes, or 
module eigengenes in the multi-set format (see \code{\link{checkSets}}). The multi-set format is a vector of
lists, one per set. Each set must contain a component \code{data} whose
rows correspond to samples and columns to eigengenes. }

  \item{setLabels}{ A vector of character strings that label sets in \code{multiME}. }

  \item{letterSubPlots}{ logical: should subplots be lettered? }

  \item{Letters}{optional specification of a sequence of letters for lettering. Defaults to "ABCD"... }

  \item{excludeGrey}{ logical: should the grey module eigengene be excluded from the plots? }

  \item{greyLabel}{ label for the grey module. Usually either "grey" or the number 0. }

  \item{plotDendrograms}{ logical: should eigengene dendrograms be plotted? }

  \item{plotHeatmaps}{ logical: should eigengene network heatmaps be plotted? }

  \item{setMargins}{ logical: should margins be set? See \code{\link[graphics]{par}}.}

  \item{marDendro}{ a vector of length 4 giving the margin setting for dendrogram plots. See
\code{\link[graphics]{par}}. If \code{setMargins} is \code{TRUE} and \code{marDendro} is not given, the
function will provide reasonable default values. }

  \item{marHeatmap}{ a vector of length 4 giving the margin setting for heatmap plots. See 
\code{\link[graphics]{par}}. If \code{setMargins} is \code{TRUE} and \code{marDendro} is not given, the  
function will provide reasonable default values. }

  \item{colorLabels}{ logical: should module eigengene names be interpreted as color names and the colors
used to label heatmap plots and barplots? }

  \item{signed}{ logical: should eigengene networks be constructed as signed? }

  \item{heatmapColors}{ color palette for heatmaps. Defaults to \code{\link{heat.colors}} when
\code{signed} is \code{FALSE}, and to \code{\link{redWhiteGreen}} when \code{signed} is \code{TRUE}. }

  \item{plotAdjacency}{ logical: should module eigengene heatmaps plot adjacency (ranging from 0 to 1),
or correlation (ranging from -1 to 1)? }

  \item{printAdjacency}{ logical: should the numerical values be printed into the adjacency or
correlation heatmap? }

  \item{cex.adjacency}{ character expansion factor for printing of numerical values into the adjacency or 
correlation heatmap }

  \item{coloredBarplot}{ logical: should the barplot of eigengene adjacency preservation distinguish
individual contributions by color? This is possible only if \code{colorLabels} is \code{TRUE} and module
eigengene names encode valid colors. }

  \item{barplotMeans}{ logical: plot mean preservation in the barplot? This option effectively rescales
the preservation by the number of eigengenes in the network. If means are plotted, the barplot is not
colored. }

  \item{barplotErrors}{ logical: should standard errors of the mean preservation be plotted? }

  \item{plotPreservation}{ a character string specifying which type of preservation measure to plot.
Allowed values are (unique abbreviations of)  \code{"standard"}, \code{"hyperbolic"}, \code{"both"}. }

  \item{zlimPreservation}{ a vector of length 2 giving the value limits for the preservation heatmaps. }

  \item{printPreservation}{ logical: should preservation values be printed within the heatmap? }

  \item{cex.preservation}{ character expansion factor for preservation display. }

  \item{\dots}{ other graphical arguments to function \code{\link{labeledHeatmap}}. }
}
\details{

  Consensus eigengene networks consist of a fixed set of eigengenes "expressed" in several different
sets. Network connection strengths are given by eigengene correlations. This function aims to visualize
the networks as well as their similarities and differences across sets. 

  The function partitions the screen appropriately and plots eigengene dendrograms in the top row, then a
square matrix of plots: heatmap plots of eigengene networks in each set on the diagonal, heatmap plots of
pairwise preservation networks below the diagonal, and barplots of aggregate network preservation of
individual eigengenes above the diagonal. A preservation plot or barplot in the row i and column j of the
square matrix represents the preservation between sets i and j. 

  Individual eigengenes are labeled by their name in the dendrograms; in the heatmaps and barplots they
can optionally be labeled by color squares. For compatibility with other functions, the color labels are
encoded in the eigengene names by prefixing the color with two letters, such as \code{"MEturquoise"}.  

  Two types of network preservation can be plotted: the \code{"standard"} is simply the difference
between adjacencies in the two compared sets. The \code{"hyperbolic"} difference de-emphasizes the
preservation of low adjacencies. When \code{"both"} is specified, standard preservation is plotted in the
lower triangle and hyperbolic in the upper triangle of each preservation heatmap. 

If the eigengenes are labeled by color, the bars in the barplot can be split into segments representing
the contribution of each eigengene and labeled by the contribution. For example, a yellow segment in a
bar labeled by a turquoise square represents the preservation of the adjacency between the yellow and
turquoise eigengenes in the two networks compared by the barplot. 

  For large numbers of eigengenes and/or sets, it may be difficult to get a meaningful plot fit a
standard computer screen. In such cases we recommend using a device such as \code{\link{postscript}} or
\code{\link{pdf}} where the user can specify large dimensions; such plots can be conveniently viewed in
standard pdf or postscript viewers.
}
\value{
  None.
}
\references{ 

  For theory and applications of consensus eigengene networks, see

  Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between co-expression
  modules. BMC Systems Biology 2007, 1:54
}

\author{ Peter Langfelder }

\seealso{ 

   \code{\link{labeledHeatmap}}, \code{\link{labeledBarplot}} for annotated heatmaps and barplots;

   \code{\link[stats]{hclust}} for hierarchical clustering and dendrogram plots

}
\keyword{ hplot }
