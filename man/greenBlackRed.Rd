\name{greenBlackRed}
\alias{greenBlackRed}
\title{ Green-black-red color sequence }
\description{
  Generate a green-black-red color sequence of a given length.
}
\usage{
greenBlackRed(n, gamma = 1)
}
\arguments{
  \item{n}{ number of colors to be returned }
  \item{gamma}{ color correction power }
}
\details{
  The function returns a color vector that starts with pure green, gradually turns into black and then to
red. The power \code{gamma} can be used to control the behaviour of the quarter- and three quarter-values
(between green and black, and black and red, respectively). Higher powers will make the mid-colors more
green and red, respectively.
}
\value{
  A vector of colors of length \code{n}.
}
\author{ Peter Langfelder }
\examples{
  par(mfrow = c(3, 1))
  displayColors(greenBlackRed(50));
  displayColors(greenBlackRed(50, 2));
  displayColors(greenBlackRed(50, 0.5));
}
\keyword{color}% __ONLY ONE__ keyword per line
