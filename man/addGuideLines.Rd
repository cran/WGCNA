\name{addGuideLines}
\alias{addGuideLines}

\title{ Add vertical ``guide lines'' to a dendrogram plot}
\description{
  Adds vertical ``guide lines'' to a dendrogram plot.
}
\usage{
addGuideLines(dendro, 
              all = FALSE, 
              count = 50, 
              positions = NULL, 
              col = "grey30", 
              lty = 3, 
              hang = 0)
}
\arguments{
  \item{dendro}{ The dendrogram (see \code{\link{hclust}}) to which the guide lines are to be added. }
  \item{all}{ Add a guide line to every object on the dendrogram? Useful if the number of objects is
relatively low. }
  \item{count}{ Number of guide lines to be plotted. The lines will be equidistantly spaced. }
  \item{positions}{ Horizontal positions of the added guide lines. If given, overrides \code{count}. }
  \item{col}{ Color of the guide lines }
  \item{lty}{ Line type of the guide lines. See \code{\link{par}}. }
  \item{hang}{ Fraction of the figure height that will separate top ends of guide lines and the
merge heights of the corresponding objects. }
}
\author{ Peter Langfelder }
\keyword{ hplot }
