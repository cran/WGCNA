\name{addTraitToMEs}
\alias{addTraitToMEs}
\title{ Add trait information to multi-set module eigengene structure }
\description{
   Adds trait information to multi-set module eigengene structure.
}
\usage{
addTraitToMEs(multiME, multiTraits)
}
\arguments{
  \item{multiME}{ Module eigengenes in multi-set format. A vector of lists, one list per set. Each list
must contain an element named \code{data} that is a data frame with module eigengenes. }
  \item{multiTraits}{ Microarray sample trait(s) in multi-set format. A vector of lists, one list per
set. Each list 
must contain an element named \code{data} that is a data frame in which each column corresponds to a
trait, and each row to an individual sample.}
}
\details{
  The function simply \code{cbind}'s the module eigengenes and traits for each set. The number of sets
and numbers of samples in each set must be consistent between \code{multiMEs} and \code{multiTraits}.
}
\value{
  A multi-set structure analogous to the input: a vector of lists, one list per set. Each list will
contain a component \code{data} with the merged eigengenes and traits for the corresponding set.
}
\author{ Peter Langfelder }
\seealso{ \code{\link{checkSets}}, \code{\link{moduleEigengenes}} }
\keyword{ misc }
