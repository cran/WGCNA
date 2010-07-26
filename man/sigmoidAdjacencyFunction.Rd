\name{sigmoidAdjacencyFunction}
\alias{sigmoidAdjacencyFunction}
\title{ Sigmoid-type adacency function. }
\description{
  Sigmoid-type function that converts a similarity to a weighted network adjacency.
}
\usage{
sigmoidAdjacencyFunction(ss, mu = 0.8, alpha = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{ similarity, a number between 0 and 1. Can be given as a scalar, vector or a matrix. }
  \item{mu}{ shift parameter. }
  \item{alpha}{ slope parameter. }
}
\details{
  The sigmoid adjacency function is defined as \eqn{1/(1+\exp[-\alpha(ss - \mu)])}{1/(1 + exp(-alpha *
(ss - mu)))}.
}
\value{

  Adjacencies returned in the same form as the input \code{ss}
  
}
\references{ Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
}
\author{ Steve Horvath }
\keyword{misc }% __ONLY ONE__ keyword per line
