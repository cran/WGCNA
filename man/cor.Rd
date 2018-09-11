\name{cor}
\alias{cor1}
\alias{corFast}
\alias{cor}
\title{ Fast calculations of Pearson correlation. }
\description{
  These functions implements a faster calculation of (weighted) Pearson correlation. 

  The speedup against the R's standard \code{\link[stats]{cor}} function will be substantial particularly
if the input matrix only contains a small number of missing data. If there are no missing data, or the
missing data are numerous, the speedup will be smaller.
}
\usage{
cor(x, y = NULL, 
    use = "all.obs", 
    method = c("pearson", "kendall", "spearman"),
    weights.x = NULL,
    weights.y = NULL,
    quick = 0, 
    cosine = FALSE, 
    cosineX = cosine,
    cosineY = cosine, 
    drop = FALSE,
    nThreads = 0, 
    verbose = 0, indent = 0)

corFast(x, y = NULL, 
    use = "all.obs", 
    quick = 0, nThreads = 0, 
    verbose = 0, indent = 0)

cor1(x, use = "all.obs", verbose = 0, indent = 0)

}
\arguments{
  \item{x}{ a numeric vector or a matrix. If \code{y} is null, \code{x} must be a matrix. }

  \item{y}{ a numeric vector or a matrix. If not given, correlations of columns of \code{x} will be
calculated. }

  \item{use}{ a character string specifying the handling of missing data. The fast calculations currently
support \code{"all.obs"} and \code{"pairwise.complete.obs"}; for other options, see R's standard
correlation function \code{\link[stats]{cor}}.  Abbreviations are allowed. }

  \item{method}{ a character string specifying the method to be used. Fast calculations are currently
available only for \code{"pearson"}. }

  \item{weights.x}{optional observation weights for \code{x}. A matrix of the same dimensions as \code{x}, 
containing non-negative weights. Only used in fast calculations: \code{methods} must be \code{"pearson"} and
\code{use} must be one of \code{"all.obs", "pairwise.complete.obs"}.}

  \item{weights.y}{optional observation weights for \code{y}. A matrix of the same dimensions as \code{y}, 
containing non-negative weights. Only used in fast calculations: \code{methods} must be \code{"pearson"} and 
\code{use} must be one of \code{"all.obs", "pairwise.complete.obs"}.}

  \item{quick}{ real number between 0 and 1 that controls the precision of handling of missing data in the
calculation of correlations. See details. }

  \item{cosine}{ logical: calculate cosine correlation? Only valid for \code{method="pearson"}.
Cosine correlation is similar to Pearson correlation but the mean subtraction is not performed. The result
is the cosine of the angle(s) between (the columns of) \code{x} and \code{y}. }

  \item{cosineX}{ logical: use the cosine calculation for \code{x}? This setting does not affect \code{y}
and can be used to give a hybrid cosine-standard correlation. }

  \item{cosineY}{ logical: use the cosine calculation for \code{y}? This setting does not affect \code{x}
and can be used to give a hybrid cosine-standard correlation. }

  \item{drop}{logical: should the result be turned into a vector if it is effectively one-dimensional? }

  \item{nThreads}{ non-negative integer specifying the number of parallel threads to be used by certain
parts of correlation calculations. 
This option only has an effect on systems on which a POSIX thread
library is available (which currently includes Linux and Mac OSX, but excludes Windows).
If zero, the number of online processors will be used if it can be determined dynamically, otherwise
correlation calculations will use 2 threads. 
Note that this option does not affect what is usually the most expensive part of the
calculation, namely the matrix multiplication. The matrix multiplication is carried out by BLAS routines provided by R;
these can be sped up by installing a fast BLAS and making R use it. }

  \item{verbose}{ Controls the level of verbosity. Values above zero will cause a small amount of
diagnostic messages to be printed. }

  \item{indent}{ Indentation of printed diagnostic messages. Each unit above zero adds two spaces.}
}

\details{
  The fast calculations are currently implemented only for \code{method="pearson"} and \code{use} either
\code{"all.obs"} or \code{"pairwise.complete.obs"}. 
  The \code{corFast} function is a wrapper that calls the function \code{cor}. If the combination of
\code{method} and \code{use} is implemented by the fast calculations, the fast code is executed; 
otherwise, R's own correlation \code{\link[stats]{cor}} is executed.



  The argument \code{quick} specifies the precision of handling of missing data. Zero will cause all
calculations to be executed precisely, which may be significantly slower than calculations without
missing data. Progressively higher values will speed up the
calculations but introduce progressively larger errors. Without missing data, all column means and
variances can be pre-calculated before the covariances are calculated. When missing data are present,
exact calculations require the column means and variances to be calculated for each covariance. The
approximate calculation uses the pre-calculated mean and variance and simply ignores missing data in the
covariance calculation. If the number of missing data is high, the pre-calculated means and variances may
be very different from the actual ones, thus potentially introducing large errors. 
The \code{quick} value times the
number of rows specifies the maximum difference in the 
number of missing entries for mean and variance calculations on the one hand and covariance on the other   
hand that will be tolerated before a recalculation is triggered. The hope is that if only a few missing
data are treated approximately, the error introduced will be small but the potential speedup can be
significant. 

}

\value{
  The matrix of the Pearson correlations of the columns of \code{x} with columns of \code{y} if \code{y}
is given, and the correlations of the columns of \code{x} if \code{y} is not given. 
}
\author{ Peter Langfelder }

\references{
  Peter Langfelder, Steve Horvath (2012)
  Fast R Functions for Robust Correlations and Hierarchical Clustering.
  Journal of Statistical Software, 46(11), 1-17.
  \url{http://www.jstatsoft.org/v46/i11/}
}

\note{ 
  The implementation uses the BLAS library matrix multiplication function for the most expensive step of
the calculation. Using a tuned, architecture-specific BLAS may significantly improve the performance of
this function.

The values returned by the corFast function may differ from the values returned by R's function
\code{\link[stats]{cor}} by rounding errors on the order of 1e-15. 

}
\seealso{ R's standard Pearson correlation function \code{\link{cor}}. }
\examples{

## Test the speedup compared to standard function cor

# Generate a random matrix with 200 rows and 1000 columns

set.seed(10)
nrow = 100;
ncol = 500;
data = matrix(rnorm(nrow*ncol), nrow, ncol);

## First test: no missing data

system.time( {corStd = stats::cor(data)} );

system.time( {corFast = cor(data)} );

all.equal(corStd, corFast)

# Here R's standard correlation performs very well.

# We now add a few missing entries.

data[sample(nrow, 10), 1] = NA;

# And test the correlations again...

system.time( {corStd = stats::cor(data, use ='p')} );

system.time( {corFast = cor(data, use = 'p')} );

all.equal(corStd, corFast)

# Here the R's standard correlation slows down considerably
# while corFast still retains it speed. Choosing
# higher ncol above will make the difference more pronounced.

}

\keyword{ misc }
