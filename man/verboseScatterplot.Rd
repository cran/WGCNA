\name{verboseScatterplot}
\alias{verboseScatterplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Scatterplot annotated by regression line and p-value}
\description{
  Produce a scatterplot annotated by the correlation, p-value, and regression line. 
}
\usage{
verboseScatterplot(x, y, 
                   sample = NULL,
                   corFnc = "cor", corOptions = "use = 'p'", 
                   main = "", xlab = NA, ylab = NA, 
                   cex = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
                   abline = FALSE, abline.color = 1, abline.lty = 1,
                   corLabel = corFnc, 
                   displayAsZero = 1e-5,
                   col = 1, bg = 0, pch = 1,
                   lmFnc = lm,
                   plotPriority = NULL,
                   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{  numerical vector to be plotted along the x axis. }

  \item{y}{  numerical vector to be plotted along the y axis. }

  \item{sample}{ determines whether \code{x} and \code{y} should be sampled for plotting, useful 
to keep the plot manageable when \code{x} and \code{y} are large vectors.
The default \code{NULL} value implies no sampling. A single numeric value will be interpreted as the
number of points to sample randomly. If a vector is given, it will be interpreted as the indices of the
entries in \code{x} and \code{y} that should be plotted. In either case, the correlation and p value will
be determined from the full vectors \code{x} and \code{y}.}

  \item{corFnc}{ character string giving the correlation function to annotate the plot. }

  \item{corOptions}{ character string giving further options to the correlation function. }

  \item{main}{ main title for the plot.}

  \item{xlab}{ label for the x-axis. }

  \item{ylab}{ label for the y-axis. }

  \item{cex}{ character expansion factor for plot annotations, recycled as necessary. }

  \item{cex.axis}{ character expansion factor for axis annotations. }

  \item{cex.lab}{ character expansion factor for axis labels. }

  \item{cex.main}{ character expansion factor for the main title. }

  \item{abline}{ logical: should the linear regression fit line be plotted? }

  \item{abline.color}{ color specification for the fit line.}

  \item{abline.lty}{ line type for the fit line.}

  \item{corLabel}{ character string to be used as the label for the correlation value printed in the main
title. }

  \item{displayAsZero}{ Correlations whose absolute value is smaller than this number will be displayed as
zero. This can result in a more intuitive display (for example, cor=0 instead of cor=2.6e-17).}

  \item{col}{color of the plotted symbols. Recycled as necessary. }

  \item{bg}{fill color of the plotted symbols (used for certain symbols). Recycled as necessary. }

  \item{pch}{Integer code for plotted symbols (see \code{link{plot.default}}). Recycled as necessary. }
  \item{lmFnc}{linear model fit function. Used to calculate the linear model fit line if \code{'abline'} is
               \code{TRUE}. For example, robust linear models are implemented in the 
               function \code{\link[MASS]{rlm}}. }

  \item{plotPriority}{Optional numeric vector of same length as \code{x}. Points with higher plot priority
will be plotted later, making them more visible if points overlap.}

  \item{\dots}{ other arguments to the function \code{\link{plot}}. }
}
\details{
  Irrespective of the specified correlation function, the p-value is always calculated for pearson
correlation. 
}
\value{
  If \code{sample} above is given, the indices of the plotted points are returned invisibly.
}
\author{ Steve Horvath and Peter Langfelder }
\seealso{ \code{\link{plot.default}} for standard scatterplots }

\keyword{hplot}% __ONLY ONE__ keyword per line
