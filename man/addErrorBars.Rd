\name{addErrorBars}
\alias{addErrorBars}
\title{ Add error bars to a barplot. }
\description{
  This function adds error bars to an existing barplot.
}
\usage{
addErrorBars(means, errors, two.side = FALSE)
}
\arguments{
  \item{means}{ vector of means plotted in the barplot }
  \item{errors}{ vector of standard errors (signle positive values) to be plotted. }
  \item{two.side}{ should the error bars be two-sided? }
}
\value{
  None.
}
\author{ Steve Horvath and Peter Langfelder }
\keyword{hplot}
