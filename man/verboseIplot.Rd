\name{verboseIplot}
\alias{verboseIplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatterplot with density 
}
\description{
Produce a scatterplot that shows density with color and is annotated by the correlation, MSE, and regression line.
}
\usage{
verboseIplot(
             x, y, 
             xlim = NA, ylim = NA, 
             nBinsX = 150, nBinsY = 150, 
             ztransf = function(x) {x}, gamma = 1, 
             sample = NULL, corFnc = "cor", corOptions = "use = 'p'", 
             main = "", xlab = NA, ylab = NA, cex = 1, 
             cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
             abline = FALSE, abline.color = 1, abline.lty = 1, 
             corLabel = corFnc, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numerical vector to be plotted along the x axis.
}
  \item{y}{
numerical vector to be plotted along the y axis.
}
  \item{xlim}{
define the range in x axis
}
  \item{ylim}{
define the range in y axis
}
  \item{nBinsX}{
number of bins along the x axis
}
  \item{nBinsY}{
number of bins along the y axis
}
  \item{ztransf}{
Function to transform the number of counts per pixel, which will be mapped by the function in colramp to well defined colors. The user has to make sure that the transformed density lies in the range [0,zmax], where zmax is any positive number (>=2). 
}
  \item{gamma}{
color correction power 
}
  \item{sample}{
  either a number of points to be sampled or a vector of indices input \code{x} and \code{y} 
for points to be plotted.  Useful when the input vectors are large and plotting all points is not practical.
}
  \item{corFnc}{
character string giving the correlation function to annotate the plot.
}
  \item{corOptions}{
character string giving further options to the correlation function.
}
  \item{main}{
main title for the plot.
}
  \item{xlab}{
label for the x-axis.
}
  \item{ylab}{
label for the y-axis.
}
  \item{cex}{
character expansion factor for plot annotations.
}
  \item{cex.axis}{
character expansion factor for axis annotations.
}
  \item{cex.lab}{
character expansion factor for axis labels.
}
  \item{cex.main}{
character expansion factor for the main title.
}
  \item{abline}{
logical: should the linear regression fit line be plotted?
}
  \item{abline.color}{
color specification for the fit line.
}
  \item{abline.lty}{
line type for the fit line.
}
  \item{corLabel}{
character string to be used as the label for the correlation value printed in the main title.
}
  \item{\dots}{
other arguments to the function plot.
}
}
\details{
Irrespective of the specified correlation function, the MSE is always calculated based on the residuals of a
linear model.
}
\value{
If sample above is given, the indices of the plotted points are returned invisibly.

}
\author{
Chaochao Cai, Steve Horvath

}
\note{
This funtion is based on verboseScatterplot (Steve Horvath and Peter Langfelder), iplot (Andreas Ruckstuhl, Rene Locher) and greenWhiteRed(Peter Langfelder )
}

\seealso{
\link{image} for more parameters
}
\keyword{graphics }
