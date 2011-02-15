\name{na}

\alias{log.na}
\alias{sum.na}
\alias{mean.na}
\alias{var.na}
\alias{cor.na}
\alias{quantile.na}
\alias{length.na}
\alias{order.na}
\alias{scale.na}
\alias{prod.na}

\title{Basic Statistical Functions for Handling Missing Values}

\description{
Basic statistical functions for handling missing values or NA. \cr 
In \code{log.na}, \code{sum.na}, \code{mean.na} and \code{var.na},
\code{quantile.na}, \code{length.na}, missing values are omitted
from the calculation. \cr 
The function \code{cor.na} calls \code{cor} with the argument
\code{use="pairwise.complete.obs"}. \cr 
The function \code{order.na} only handles vector arguments and not
lists.  However, it gives the option of omitting the NAs
(\code{na.last=NA}), of placing the NAs at the start of the ordered
vector (\code{na.last=F}) or at the end (\code{na.last=T}). \cr 
The function \code{scale.na} is a modified version of
\code{\link{scale}} which allows NAs in the variance calculation. If
\code{scale = T}, the function \code{f} in \code{scale.na} uses
\code{var.na} to perform the variance calculation.
The function \code{prod.na} is similar to the \code{\link{prod}}
function with \code{na.rm=TRUE}. This function returns the product of
all the values present in its arguments, omitting any missing values.
}

\author{
  Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
}

\seealso{\code{\link{log}}, \code{\link{sum}}, \code{\link{mean}},
  \code{\link{var}}, \code{\link{cor}}, \code{\link{order}},
  \code{\link{scale}}, \code{link{prod}}.}

\keyword{misc}

