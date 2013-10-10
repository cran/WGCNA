\name{simulateEigengeneNetwork}
\alias{simulateEigengeneNetwork}
\title{ Simulate eigengene network from a causal model }
\description{
  Simulates a set of eigengenes (vectors) from a given set of causal anchors and a causal matrix.
}
\usage{
simulateEigengeneNetwork(
  causeMat, 
  anchorIndex, anchorVectors, 
  noise = 1, 
  verbose = 0, indent = 0)
}
\arguments{
  \item{causeMat}{ causal matrix. The entry \code{[i,j]} is the influence (path coefficient) of
vector \code{j} on vector \code{i}. }  
  \item{anchorIndex}{ specifies the indices of the anchor vectors. }
  \item{anchorVectors}{ a matrix giving the actual anchor vectors as columns. Their number must
equal the length of \code{anchorIndex}. }
  \item{noise}{ standard deviation of the noise added to each simulated vector. }
  \item{verbose}{ level of verbosity. 0 means silent. }
  \item{indent}{ indentation for diagnostic messages. Zero means no indentation; each unit adds two
spaces. }
}
\details{
  The algorithm starts with the anchor vectors and iteratively generates the rest from the path
coefficients given in the matrix \code{causeMat}.
}
\value{
  A list with the following components:
  \item{eigengenes }{ generated eigengenes. }
  \item{causeMat }{ a copy of the input causal matrix}
  \item{levels}{ useful for debugging. A vector with one entry for each eigengene giving the number
of generations of parents of the eigengene. Anchors have level 0, their direct causal children have
level 1 etc.}
  \item{anchorIndex}{a copy of the input \code{anchorIndex}. } 
}
\author{ Peter Langfelder }
\keyword{ misc }% __ONLY ONE__ keyword per line
