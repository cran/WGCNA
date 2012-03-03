\name{TOMplot}
\alias{TOMplot}
\title{ Graphical representation of the Topological Overlap Matrix }
\description{
  Graphical representation of the Topological Overlap Matrix using a heatmap plot combined with
the corresponding hierarchical clustering dendrogram and module colors. 
}
\usage{
TOMplot(
   dissim, 
   dendro, 
   Colors = NULL, 
   ColorsLeft = Colors, 
   terrainColors = FALSE, 
   setLayout = TRUE,
   ...)
}
\arguments{
  \item{dissim}{ a matrix containing the topological overlap-based dissimilarity }

  \item{dendro}{ the corresponding hierarchical clustering dendrogram }

  \item{Colors}{ optional specification of module colors to be plotted on top }

  \item{ColorsLeft}{ optional specification of module colors on the left side. If \code{NULL},
\code{Colors} will be used.  }

  \item{terrainColors}{ logical: should terrain colors be used? }

  \item{setLayout}{ logical: should layout be set? If \code{TRUE}, standard layout for one plot will be
used. Note that this precludes multiple plots on one page. If \code{FALSE}, the user is responsible for
setting the correct layout. }

  \item{\dots}{ other graphical parameters to \code{\link{heatmap}}. }
}
\details{
  The standard \code{heatmap} function uses the \code{\link{layout}} function to set the following
layout (when \code{Colors} is given):
\preformatted{
0 0 5
0 0 2
4 1 3
}
To get a meaningful heatmap plot, user-set layout must respect this geometry. 
}
\value{
  None.
}
\author{ Steve Horvath and Peter Langfelder }
\seealso{ \code{\link{heatmap}}, the workhorse function doing the plotting. }
\keyword{ misc }
