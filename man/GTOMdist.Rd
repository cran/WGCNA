\name{GTOMdist}
\alias{GTOMdist}
\title{ Generalized Topological Overlap Measure }
\description{
  Generalized Topological Overlap Measure, taking into account interactions of higher degree.
}
\usage{
GTOMdist(adjMat, degree = 1)
}
\arguments{
  \item{adjMat}{ adjacency matrix. See details below. }
  \item{degree}{ integer specifying the maximum degree to be calculated. }
}
\value{
  Matrix of the same dimension as the input \code{adjMat}.
}
\references{ Yip A, Horvath S (2007) Gene network interconnectedness and the generalized topological
overlap measure. BMC Bioinformatics 2007, 8:22 }
\author{ Steve Horvath and Andy Yip }
\keyword{ misc }% __ONLY ONE__ keyword per line
