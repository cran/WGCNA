\name{stat.bwss}

\alias{stat.bwss}


\title{Between and Within Group Sum of Squares Calculation}

\description{This function computes the between and within group sum
of squares for each row of a matrix which may have NAs. 
}

\usage{
stat.bwss(x, cl)
}

\arguments{
 \item{x}{a matrix, NAs allowed. In the microarray context, the rows
 of \code{X} could correspond to genes and the columns to different
 hybridizations.} 
 \item{cl}{column labels, must be consecutive integers.} 
}

\value{
List containing the following components
  \item{wn}{matrix with class sizes for each row of \code{X};}
  \item{BW}{vector of BSS/WSS for each row of \code{X};}
  \item{BSS}{between group sum of squares for each row of \code{X};}
  \item{WSS}{within group sum of squares for each row of \code{X};}
  \item{TSS}{total sum of squares for each row of \code{X};}
  \item{tvar}{variance for each row of \code{X}.}
}

\references{S. Dudoit, J. Fridlyand, and T. P. Speed. Comparison of
Discrimination Methods for the Classification of Tumors Using Gene
Expression Data. June 2000. (Statistics, UC Berkeley,  Tech Report \#
576). } 

\author{
  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}  \cr
  Jane Fridlyand, \email{janef@stat.berkeley.edu}
}

\seealso{\code{\link{var.na}}, \code{\link{var}}, \code{\link{apply}}.}

\keyword{misc}
