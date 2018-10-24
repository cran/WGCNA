\name{consensusTreeInputs}
\alias{consensusTreeInputs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get all elementary inputs in a consensus tree
}
\description{
This function returns a flat vector or a structured list of elementary inputs to a given consensus tree, that is, inputs
that are not consensus trees themselves.
}
\usage{
consensusTreeInputs(consensusTree, flatten = TRUE)
}
\arguments{
  \item{consensusTree}{
A consensus tree of class \code{\link{ConsensusTree}}.
}
  \item{flatten}{
Logical; if \code{TRUE}, the function returns a flat character vector of inputs; otherwise, a list whose structure reflects 
the structure of \code{consensusTree}.
}
}
\value{
A character vector of inputs or a list of inputs whose structure reflects 
the structure of \code{consensusTree}.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{newConsensusTree}} for creating consensus trees.
}
\keyword{misc}% __ONLY ONE__ keyword per line
