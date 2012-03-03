\name{bicorAndPvalue}
\alias{bicorAndPvalue}
\title{
Calculation of biweight midcorrelations and associated p-values
}
\description{
 A faster, one-step calculation of Student correlation p-values for multiple biweight 
midcorrelations, properly taking
into account the actual number of observations. 
}
\usage{
bicorAndPvalue(x, y = NULL, 
             use = "pairwise.complete.obs", 
             alternative = c("two.sided", "less", "greater"),
             ...)
}
\arguments{
  \item{x}{
a vector or a matrix
}
  \item{y}{
a vector or a matrix. If \code{NULL}, the correlation of columns of \code{x} will be calculated.
}
  \item{use}{
determines handling of missing data. See \code{\link{bicor}} for details.
}
  \item{alternative}{
specifies the alternative hypothesis and must be (a unique abbreviation of) one of
\code{"two.sided"}, \code{"greater"} or \code{"less"}.
the initial letter.  \code{"greater"} corresponds to positive
association, \code{"less"} to negative association.
}
  \item{\dots}{
other arguments to the function \code{\link{bicor}}. 
}
}
\details{
The function calculates the biweight midcorrelations of a matrix or of two matrices 
and the corresponding Student p-values.
The output is not as full-featured as \code{\link{cor.test}}, but can work with matrices as input.
}
\value{
A list with the following components, each a marix:
 \item{bicor}{the calculated correlations}
 \item{p}{the Student p-values corresponding to the calculated correlations}
 \item{Z}{Fisher transform of the calculated correlations}
 \item{t}{Student t statistics of the calculated correlations}
 \item{nObs}{Numbers of observations for the correlation, p-values etc.}
}
\author{
Peter Langfelder and Steve Horvath 
}

\references{
  Peter Langfelder, Steve Horvath (2012)
  Fast R Functions for Robust Correlations and Hierarchical Clustering.
  Journal of Statistical Software, 46(11), 1-17.
  \url{http://www.jstatsoft.org/v46/i11/}
}

\seealso{
\code{\link{bicor}} for calculation of correlations only;

\code{\link{cor.test}} for another function for significance test of correlations
}
\examples{
# generate random data with non-zero correlation
set.seed(1);
a = rnorm(100);
b = rnorm(100) + a;
x = cbind(a, b);
# Call the function and display all results
bicorAndPvalue(x)
# Set some components to NA
x[c(1:4), 1] = NA
corAndPvalue(x)
# Note that changed number of observations.
}
\keyword{ stats }
