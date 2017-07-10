\name{plotMultiHist}
\alias{plotMultiHist}
\title{
Plot multiple histograms in a single plot 
}
\description{
This function plots density or cumulative distribution function of multiple histograms in a single plot, using
lines.
}
\usage{
plotMultiHist(
   data, 
   nBreaks = 100, 
   col = 1:length(data), 
   scaleBy = c("area", "max", "none"), 
   cumulative = FALSE, 
   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A list in which each component corresponds to a separate histogram and is a vector of values to be shown in
each histogram.
}
  \item{nBreaks}{
Number of breaks in the combined plot.
}
  \item{col}{
Color of the lines. Should be a vector of the same length as \code{data}.
}
  \item{scaleBy}{
Method to make the different histograms comparable. The counts are scaled such that either the total area or
the maximum are the same for all histograms, or the histograms are shown without scaling.
}
  \item{cumulative}{
Logical: should the cumulative distribution be shown instead of the density?
}
  \item{\dots}{
Other graphical arguments.
}
}
\value{
  Invisibly, 
\item{x}{A list with one component per histogram (component of \code{data}), giving the bin midpoints}
\item{y}{A list with one component per histogram (component of \code{data}), giving the scaled bin counts}
}
\author{
Peter Langfelder
}
\note{
This function is still experimental and behavior may change in the future.
}

\seealso{
\code{\link{hist}}
}
\examples{
data = list(rnorm(1000), rnorm(10000) + 2);
plotMultiHist(data, xlab = "value", ylab = "scaled density")
}
\keyword{misc}% __ONLY ONE__ keyword per line
