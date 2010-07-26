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
                   abline = FALSE, corLabel = corFnc, ...)
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

  \item{cex}{ character expansion factor for plot annotations. }

  \item{cex.axis}{ character expansion factor for axis annotations. }

  \item{cex.lab}{ character expansion factor for axis labels. }

  \item{cex.main}{ character expansion factor for the main title. }

  \item{abline}{ logical: should the linear regression fit line be plotted? }

  \item{corLabel}{ character string to be used as the label for the correlation value printed in the main
title. }

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
