\name{randIndex}
\alias{randIndex}
\title{ Rand index of two partitions}
\description{
  Computes the Rand index, a measure of the similarity between two clusterings.
}
\usage{
randIndex(tab, adjust = TRUE)
}
\arguments{
  \item{tab}{ a matrix giving the cross-tabulation table of two clusterings. }
  \item{adjust}{logical: should the "adjusted" version be computed? }
}
\value{
  the Rand index of the input table.
}
\references{ W. M. Rand (1971). "Objective criteria for the evaluation of clustering methods". Journal of the
American Statistical Association 66: 846-850}
\author{ Steve Horvath}
\keyword{ misc }
