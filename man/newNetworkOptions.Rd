\name{newNetworkOptions}
\alias{newNetworkOptions}
\alias{NetworkOptions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a list of network construction arguments (options).
}
\description{
This function creates a reusable list of network calculation arguments/options.
}
\usage{
newNetworkOptions(
    correlationOptions = newCorrelationOptions(),

    # Adjacency options
    replaceMissingAdjacencies = TRUE,
    power = 6,
    networkType = c("signed hybrid", "signed", "unsigned"),
    checkPower = TRUE,

    # Topological overlap options
    TOMType = c("signed", "signed Nowick", "unsigned", "none",
                "signed 2", "signed Nowick 2", "unsigned 2"),
    TOMDenom = c("mean", "min"),
    suppressTOMForZeroAdjacencies = FALSE,
    suppressNegativeTOM = FALSE,

    # Internal behavior options
    useInternalMatrixAlgebra = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{correlationOptions}{
A list of correlation options. See \code{\link{newCorrelationOptions}}.
}
  \item{replaceMissingAdjacencies}{Logical: should missing adjacencies be replaced by zero? 
}
  \item{power}{ Soft-thresholding power for network construction. 
}
  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. 
}
  \item{checkPower}{Logicel: should the power be checked for sanity? 
}
  \item{TOMType}{One of \code{"none"}, \code{"unsigned"}, \code{"signed"}, \code{"signed Nowick"},
\code{"unsigned 2"}, \code{"signed 2"} and \code{"signed Nowick 2"}. If \code{"none"}, adjacency
will be used for clustering. See \code{\link{TOMsimilarityFromExpr}} for details.}

  \item{TOMDenom}{Character string specifying the TOM variant to be used. Recognized values are
\code{"min"} giving the standard TOM described in Zhang and Horvath (2005), and \code{"mean"} in which
the \code{min} function in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
better results but at this time should be considered experimental.
}
  \item{suppressTOMForZeroAdjacencies}{logical: for those components that have zero adjacency, should TOM be
set to zero as well?}

 \item{suppressNegativeTOM}{Logical: should the result be set to zero when negative? Negative TOM values can occur when
\code{TOMType} is \code{"signed Nowick"}.}


  \item{useInternalMatrixAlgebra}{logical: should internal implementation of matrix multiplication be used
instead of R-provided BLAS? The internal implementation is slow and this option should only be used if one
suspects a bug in R-provided BLAS.}

}
\value{
A list of class \code{NetworkOptions}.
}
\author{
Peter Langfelder
}

\seealso{
code{\link{newCorrelationOptions}}
}
\keyword{misc}
