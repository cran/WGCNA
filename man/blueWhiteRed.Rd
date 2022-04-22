\name{blueWhiteRed}
\alias{blueWhiteRed}
\title{ Blue-white-red color sequence }
\description{
  Generate a blue-white-red color sequence of a given length.
}
\usage{
blueWhiteRed(
    n, 
    gamma = 1, 
    endSaturation = 1,
    blueEnd = c(0.05 + (1-endSaturation) * 0.45 , 0.55 + (1-endSaturation) * 0.25, 1.00),
    redEnd = c(1.0, 0.2 + (1-endSaturation) * 0.6, 0.6*(1-endSaturation)),
    middle = c(1,1,1))
}
\arguments{
  \item{n}{ number of colors to be returned. }
  \item{gamma}{ color change power. }
  \item{endSaturation}{ a number between 0 and 1 giving the saturation of the colors that will represent the
ends of the scale. Lower numbers mean less saturation (lighter colors).}
  \item{blueEnd}{vector of length 3 giving the RGB relative values (between 0 and 1) for the blue or negative end color.}
  \item{redEnd}{vector of length 3 giving the RGB relative values (between 0 and 1) for the red or positive end color.}
  \item{middle}{vector of length 3 giving the RGB relative values (between 0 and 1) for the middle of the scale.}
}
\details{
  The function returns a color vector that starts with blue, gradually turns into white and then to
red. The power \code{gamma} can be used to control the behaviour of the quarter- and three quarter-values
(between blue and white, and white and red, respectively). Higher powers will make the mid-colors more
white, while lower powers will make the colors more saturated, respectively.
}
\value{
  A vector of colors of length \code{n}.
}
\author{ Peter Langfelder }
\seealso{
\code{\link{numbers2colors}} for a function that produces a color representation for continuous numbers.
}

\examples{
  par(mfrow = c(3, 1))
  displayColors(blueWhiteRed(50));
  title("gamma = 1")
  displayColors(blueWhiteRed(50, 3));
  title("gamma = 3")
  displayColors(blueWhiteRed(50, 0.5));
  title("gamma = 0.5")
}
\keyword{color}% __ONLY ONE__ keyword per line
