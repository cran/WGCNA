\name{stat.diag.da}

\alias{stat.diag.da}

\title{Diagonal Discriminant Analysis}

\description{
This function implements a simple Gaussian maximum likelihood
discriminant rule, for diagonal class covariance matrices.  
}

\usage{
stat.diag.da(ls, cll, ts, pool=1)
}

\arguments{
 \item{ls}{learning set data matrix, with rows corresponding to
cases (i.e., mRNA samples) and columns to predictor variables 
(i.e., genes).} 
 \item{cll}{class labels for learning set, must be consecutive integers.}
 \item{ts}{test set data matrix, with rows corresponding to cases
and columns to predictor variables.} 
 \item{pool}{logical flag. If \code{pool=1}, the covariance matrices
are assumed to be constant across classes and the discriminant rule
is linear in the data. If \code{pool=0}, the covariance matrices may
vary across classes and the discriminant rule is quadratic in the
data.} 
} 

\value{
 List containing the following components 
  \item{pred}{vector of class predictions for the test set.}
}

\references{
S. Dudoit, J. Fridlyand, and T. P. Speed. Comparison of
Discrimination Methods for the Classification of Tumors Using Gene 
Expression Data. June 2000. (Statistics, UC Berkeley, Tech Report \#576).}

\author{
  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
  Jane Fridlyand, \email{janef@stat.berkeley.edu}
}

%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\keyword{misc}


