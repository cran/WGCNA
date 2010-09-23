\name{pquantile}
\alias{pquantile}
\alias{pmedian}
\alias{pmean}
\title{ Parallel quantile, median, mean }
\description{
Calculation of ``parallel'' quantiles, medians, and means, across given arguments.
}
\usage{
pquantile(prob, ...)
pmedian(...)
pmean(...)
}
\arguments{
  \item{prob}{
    A number or vector of probabilities at which to calculate the quantile. See \code{\link{quantile}}. 
}
  \item{\dots}{
    Numeric arguments. All arguments must have the same dimensions. See details.
}
}
\details{
    Given the argumens, say x,y,z, of equal dimensions, the \code{pquantile} 
   calculates and returns the quantile of the first components of x,y,z, then the second components, etc.
Similarly, \code{pmedian} and \code{pmean} calculate the median and mean, respectively.
}
\value{
    A vector or array containing the quantiles, medians, or
means. The dimensions are determined by the dimensions of the input arguments and whether the \code{prob}
input is scalar or a vector.  
If any of the input variables have \code{dimnames}, the first non-NULL dimnames are copied into the output.
}
\author{
 Peter Langfelder and Steve Horvath
}

\seealso{
  \code{\link{pmin}} and \code{\link{pmax}} for analogous functions for minimum and maximum,

  \code{\link{quantile}}, \code{\link{median}}, \code{\link{mean}} for the underlying statistics.
}
\examples{

# Generate 2 simple matrices
a = matrix(c(1:12), 3, 4);
b = a+ 1;
c = a + 2;

# Set the colnames on matrix a

colnames(a) = spaste("col_", c(1:4));

# Example use

pquantile(prob = 0.5, a, b, c)
pquantile(prob = c(0, 0.5, 1), a,b, c);

pmean(a,b,c)
pmedian(a,b,c)

}
\keyword{ misc }
