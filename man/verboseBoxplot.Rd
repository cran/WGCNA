\name{verboseBoxplot}
\alias{verboseBoxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Boxplot annotated by a Kruskal-Wallis p-value}
\description{
  Plot a boxplot annotated by the Kruskal-Wallis p-value. Uses the function \code{\link[graphics]{boxplot}}
for the actual drawing.
}
\usage{
verboseBoxplot(x, g, main = "", xlab = NA, ylab = NA, 
               cex = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
               notch = TRUE, varwidth = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ numerical vector of data whose group means are to be plotted }
  \item{g}{ a factor or a an object coercible to a factor giving the groups that will go into each box.}
  \item{main}{ main title for the plot.}
  \item{xlab}{ label for the x-axis. }
  \item{ylab}{ label for the y-axis. }
  \item{cex}{ character expansion factor for plot annotations. }
  \item{cex.axis}{ character expansion factor for axis annotations. }
  \item{cex.lab}{ character expansion factor for axis labels. }
  \item{cex.main}{ character expansion factor for the main title. }
  \item{notch}{logical: should the notches be drawn? See \code{\link[graphics]{boxplot}} and
\code{\link{boxplot.stats}} for details. }
  \item{varwidth}{logical: if \code{TRUE}, the boxes are drawn with widths
          proportional to the square-roots of the number of
          observations in the groups.}
  \item{\dots}{ other arguments to the function \code{\link{boxplot}}. Of note is the argument \code{las}
that specifies label orientation. Value \code{las=1} will result in horizontal labels (the default), while
\code{las=2} will result in vertical labels, useful when the labels are long.}
}
\value{
  Returns the value returned by the function \code{\link{boxplot}}.  
}
\author{ Steve Horvath }
\seealso{ \code{\link{boxplot}} }
\keyword{ misc }% __ONLY ONE__ keyword per line
