\name{greenWhiteRed}
\alias{greenWhiteRed}
\title{ Green-white-red color sequence }
\description{
  Generate a green-white-red color sequence of a given length.
}
\usage{
greenWhiteRed(n, gamma = 1)
}
\arguments{
  \item{n}{ number of colors to be returned }
  \item{gamma}{ color correction power }
}
\details{
  The function returns a color vector that starts with pure green, gradually turns into white and then to
red. The power \code{gamma} can be used to control the behaviour of the quarter- and three quarter-values
(between green and white, and white and red, respectively). Higher powers will make the mid-colors more
white, while lower powers will make the colors more saturated, respectively.
}
\value{
  A vector of colors of length \code{n}.
}
\author{ Peter Langfelder }
\examples{
  par(mfrow = c(3, 1))
  displayColors(greenWhiteRed(50));
  displayColors(greenWhiteRed(50, 3));
  displayColors(greenWhiteRed(50, 0.5));
}
\keyword{color}% __ONLY ONE__ keyword per line
