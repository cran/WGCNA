\name{addGrid}
\alias{addGrid}

\title{ Add grid lines to an existing plot. }
\description{
  This function adds horizontal and/or vertical grid lines to an existing plot. The grid lines are
aligned with tick marks.
}
\usage{
addGrid(linesPerTick = NULL, horiz = TRUE, vert = FALSE, col = "grey30", lty = 3)
}
\arguments{
  \item{linesPerTick}{ Number of lines between successive tick marks (including the line on the
tickmarks themselves) }
  \item{horiz}{ Draw horizontal grid lines? }
  \item{vert}{ Draw vertical tick lines? }
  \item{col}{ Specifies color of the grid lines }
  \item{lty}{ Specifies line type of grid lines. See \code{\link{par}}. }
}
\details{
  If \code{linesPerTick} is not specified, it is set to 5 if number of tick s is 5 or less, and it
is set to 2 if number of ticks is greater than 5.
}
\author{ Peter Langfelder }
\note{ The function does not work whenever logarithmic scales are in use. }
\examples{
  plot(c(1:10), c(1:10))
  addGrid();
}
\keyword{hplot}% __ONLY ONE__ keyword per line
