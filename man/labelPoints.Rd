\name{labelPoints}
\Rdversion{1.1}
\alias{labelPoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Label scatterplot points 
}
\description{
  Given scatterplot point coordinates, the function tries to place labels near the points such that the
labels overlap as little as possible. User beware: the algorithm implemented here is quite primitive and
while it will help in many cases, it is by no means perfect. Consider this function experimental. We hope to
improve the algorithm in the future to make it useful in a broader range of situations.
}
\usage{
labelPoints(
   x, y, labels, 
   cex = 0.7, offs = 0.01, xpd = TRUE, 
   jiggle = 0, protectEdges = TRUE, 
   doPlot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a vector of x coordinates of the points
}
  \item{y}{
  a vector of y coordinates of the points
}
  \item{labels}{
  labels to be placed next to the points
}
  \item{cex}{
  character expansion factor for the labels
}
  \item{offs}{
  offset of the labels from the plotted coordinates in inches
}
  \item{xpd}{
  logical: controls truncating labels to fit within the plotting region. See \code{\link{par}}.
}
  \item{jiggle}{
  amount of random noise to be added to the coordinates. This may be useful if the scatterplot is too
regular (such as all points on one straight line). 
}

  \item{protectEdges}{ logical: should labels be shifted inside the (actual or virtual) frame of the plot? }

  \item{doPlot}{logical: should the labels be actually added to the plot? Value \code{FALSE} may be useful if
the user would like to simply compute the best label positions the function can come up with.}

  \item{\dots}{
  other arguments to function \code{\link{text}}. 
}
}
\details{
  The algorithm basically works by finding the direction of most surrounding points, and attempting to place
the label in the opposite direction. There are (not uncommon) situations in which this placement is
suboptimal; the author promises to further develop the function sometime in the future.

Note that this function does not plot the actual scatterplot; only the labels are plotted. Plotting the
scatterplot is the responsibility of the user. 

The argument \code{offs} needs to be carefully tuned to the size of the plotted symbols. Sorry, no automation
here yet.

The argument \code{protectEdges} can be used to shift labels that would otherwise extend beyond the plot to
within the plot. Sometimes this may cause some overlapping with other points or labels; use with care.

}
\value{
  Invisibly, a data frame with 3 columns, giving the x and y positions of the labels, and the labels themselves.
}
\author{
 Peter Langfelder
}

\seealso{
  \code{\link{plot.default}}, \code{\link{text}}
}
\examples{
# generate some random points
   set.seed(11);
   n = 20;
   x = runif(n);
   y = runif(n);

# Create a basic scatterplot
   col = standardColors(n);
   plot(x,y, pch = 21, col =1, bg = col, cex = 2.6, 
        xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.0));
   labelPoints(x, y, paste("Pt", c(1:n), sep=""), offs = 0.10, cex = 1);

# label points using longer text labels. Note the positioning is not perfect, but close enough.

   plot(x,y, pch = 21, col =1, bg = col, cex = 2.6, 
        xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.0));
   labelPoints(x, y, col, offs = 0.10, cex = 0.8);
}
\keyword{ plot }% __ONLY ONE__ keyword per line
