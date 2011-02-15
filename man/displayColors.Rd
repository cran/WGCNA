\name{displayColors}
\alias{displayColors}
\title{ Show colors used to label modules }
\description{
  The function plots a barplot using colors that label modules.
}
\usage{
displayColors(colors = NULL)
}
\arguments{
  \item{colors}{colors to be displayed. Defaults to all colors available for module labeling. }
}
\details{
  To see the first \code{n} colors, use argument \code{colors = standardColors(n)}.
}
\value{
  None.
}
\author{ Peter Langfelder }
\seealso{ \code{\link{standardColors}} }
\examples{
  displayColors(standardColors(10))
}
\keyword{ misc }
