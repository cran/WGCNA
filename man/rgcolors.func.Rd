\name{rgcolors.func}

\alias{rgcolors.func}

\title{Red and Green Color Specification}

\description{
This function creates a vector of n ``contiguous'' colors,
corresponding to n intensities (between 0 and 1) of the red, green
and blue primaries, with the blue intensities set to zero. The
values returned by \code{rgcolors.func} can be used with a
\code{col=} specification in graphics functions or in
\code{\link{par}}.  
}

\usage{
rgcolors.func(n=50)
}

\arguments{
 \item{n}{the number of colors (>= 1) to be used in the red and
 green palette. } 

}
\value{a character vector of color names. Colors are specified
directly in terms of their RGB components with a string of the form
"\#RRGGBB", where each of the pairs RR, GG, BB consist of two
hexadecimal digits giving a value in the range 00 to FF. 
 }


\author{
  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
  Jane Fridlyand, \email{janef@stat.berkeley.edu}
}

\seealso{\code{\link{plotCor}}, \code{\link{plotMat}},
\code{\link{colors}}, \code{\link{rgb}}, \code{\link{image}}.} 

\examples{
rgcolors.func(n=5)
## The following vector is returned:
## "#00FF00" "#40BF00" "#808000" "#BF4000" "#FF0000"
}

\keyword{color}

