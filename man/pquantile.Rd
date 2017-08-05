\name{pquantile}
\alias{pquantile}
\alias{pquantile.fromList}
\alias{pmedian}
\alias{pmean}
\alias{pmean.fromList}
\alias{pminWhich.fromList}
\title{ Parallel quantile, median, mean }
\description{
Calculation of ``parallel'' quantiles, minima, maxima, medians, and means, across given arguments or across lists
}
\usage{
pquantile(prob, ...)
pquantile.fromList(dataList, prob)
pmedian(...)
pmean(..., weights = NULL)
pmean.fromList(dataList, weights = NULL)
pminWhich.fromList(dataList)
}
\arguments{
  \item{prob}{
    A single probability at which to calculate the quantile. See \code{\link{quantile}}. 
}
  \item{dataList}{A list of numeric vectors or arrays, all of the same length and dimensions, over which to
calculate ``parallel'' quantiles.}
  \item{weights}{Optional vector of the same length as \code{dataList}, giving the weights to be used in the
weighted mean. If not given, unit weights will be used.}
  \item{\dots}{
    Numeric arguments. All arguments must have the same dimensions. See details.
}
}
\details{
    Given numeric arguments, say x,y,z, of equal dimensions (and length), the \code{pquantile} 
   calculates and returns the quantile of the first components of x,y,z, then the second components, etc.
Similarly, \code{pmedian} and \code{pmean} calculate the median and mean, respectively.
The funtion \code{pquantile.fromList} is identical to \code{pquantile} except that the argument
\code{dataList} replaces the ... in holding the numeric vectors over which to calculate the quantiles.
}
\value{
    \item{pquantile, pquantile.fromList}{A vector or array containing quantiles.}
    \item{pmean, pmean.fromList}{A vector or array containing means. }
    \item{pmedian}{A vector or array containing  medians.}
    \item{pminWhich.fromList}{A list with two components: \code{min} gives the minima, \code{which} gives the
indices of the elements that are the minima.}
 
Dimensions are copied from dimensions of the input arguments.
If any of the input variables have \code{dimnames}, the first non-NULL dimnames are copied into the output.
}
\author{
 Peter Langfelder and Steve Horvath
}
\seealso{
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

pmean(a,b,c)
pmedian(a,b,c)

}
\keyword{ misc }
