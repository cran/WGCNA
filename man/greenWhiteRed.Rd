\name{greenWhiteRed}
\alias{greenWhiteRed}
\title{ Green-white-red color sequence }
\description{
  Generate a green-white-red color sequence of a given length.
}
\usage{
greenWhiteRed(n, gamma = 1, warn = TRUE)
}
\arguments{
  \item{n}{ number of colors to be returned }
  \item{gamma}{ color change power }
  \item{warn}{logical: should the user be warned that this function produces a palette unsuitable for
people with most common color blindness?}
}
\details{
  The function returns a color vector that starts with green, gradually turns into white and then to
red. The power \code{gamma} can be used to control the behaviour of the quarter- and three quarter-values
(between green and white, and white and red, respectively). Higher powers will make the mid-colors more
white, while lower powers will make the colors more saturated, respectively.

  Typical use of this function is to produce (via function \code{\link{numbers2colors}}) 
a color representation of numbers within a symmetric interval
around 0, for example, the interval [-1, 1]. Note though that since green and red are not distinguishable by
people with the most common type of color blindness, we recommend using the analogous palette returned by
the function \code{\link{blueWhiteRed}}.

}
\value{
  A vector of colors of length \code{n}.
}
\author{ Peter Langfelder }
\seealso{ 
\code{\link{blueWhiteRed}} for a color sequence more friendly to people with the most common type of color
blindness;

\code{\link{numbers2colors}} for a function that produces a color representation for continuous numbers.
}
\examples{
  par(mfrow = c(3, 1))
  displayColors(greenWhiteRed(50));
  title("gamma = 1")
  displayColors(greenWhiteRed(50, 3));
  title("gamma = 3")
  displayColors(greenWhiteRed(50, 0.5));
  title("gamma = 0.5")
}
\keyword{color}% __ONLY ONE__ keyword per line
