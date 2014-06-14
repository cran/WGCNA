\name{exportNetworkToVisANT}
\alias{exportNetworkToVisANT}
\title{ Export network data in format readable by VisANT}
\description{
  Exports network data in a format readable and displayable by the VisANT software.
}
\usage{
exportNetworkToVisANT(
  adjMat, 
  file = NULL, 
  weighted = TRUE, 
  threshold = 0.5, 
  maxNConnections = NULL,
  probeToGene = NULL)
}
\arguments{

  \item{adjMat}{ adjacency matrix of the network to be exported. }

  \item{file}{ character string specifying the file name of the file in which the data should be written.
If not given, no file will be created. The file is in a plain text format. }

  \item{weighted}{ logical: should the exported network by weighted? }

  \item{threshold}{ adjacency threshold for including edges in the output. }

  \item{maxNConnections}{maximum number of exported adjacency edges. This can be used as another filter on
the exported edges.}

  \item{probeToGene}{ optional specification of a conversion between probe names (that label columns and
rows of \code{adjacency}) and gene names (that should label nodes in the output). }

}
\details{
  The adjacency matrix is checked for validity. The entries can be negative, however. The adjacency
matrix is expected to also have valid \code{names} or \code{dimnames[[2]]} that represent the probe names
of the corresponding edges. 

  Whether the output is a weighted network or not, only edges whose (absolute value of) adjacency are
above \code{threshold} will be included in the output. If \code{maxNConnections} is given, at most
\code{maxNConnections} will be included in the output.

If \code{probeToGene} is given, it is expected to have two columns, the first one corresponding to the
probe names, the second to their corresponding gene names that will be used in the output. 
}

\value{
  A data frame containing the network information suitable as input to VisANT. The same data frame is
also written into a file specified by \code{file}, if given. 
}

\references{ VisANT software is available from http://visant.bu.edu/. }

\author{ Peter Langfelder }

\keyword{ misc }% __ONLY ONE__ keyword per line
