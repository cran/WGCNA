\name{BD.getData}
\alias{BD.actualFileNames}
\alias{BD.nBlocks}
\alias{BD.blockLengths}
\alias{BD.getMetaData}
\alias{BD.getData}
\alias{BD.checkAndDeleteFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Various basic operations on \code{BlockwiseData} objects.
}
\description{
These functions implement basic operations on \code{\link{BlockwiseData}} objects.
Blockwise here means
that the data is too large to be loaded or processed in one piece and is therefore split into blocks that can
be handled one by one in a divide-and-conquer manner.
}
\usage{
BD.actualFileNames(bwData)
BD.nBlocks(bwData)
BD.blockLengths(bwData)
BD.getMetaData(bwData, blocks = NULL, simplify = TRUE)
BD.getData(bwData, blocks = NULL, simplify = TRUE)
BD.checkAndDeleteFiles(bwData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bwData}{
A \code{BlockwiseData} object.
}
  \item{blocks}{
Optional vector of integers specifying the blocks on which to execute the operation.
}
  \item{simplify}{
Logical: if the \code{blocks} argument above is of length 1, should the returned list be simplified by
removing the redundant outer \code{list} structure? 
}
}
\details{
Several functions in this package use the concept of blockwise, or "divide-and-conquer", analysis. The
BlockwiseData class is meant to hold the blockwise data, or all necessary information about blockwise data
that is saved in disk files. 

}
\value{
\item{BD.actualFileNames}{returns a vector of character strings giving the file names in which the files are
saved, or \code{NULL} if the data are held in-memory.}

\item{BD.nBlocks}{returns the number of blocks in the input object.}

\item{BD.blockLengths}{returns the block lengths (results of applying \code{\link{length}} to the data in
each block).}

\item{BD.getMetaData}{returns a list with one component per block. Each component is in turn a list
containing the stored meta-data for the corresponding block. If \code{blocks} is of length 1 and
\code{simplify} is \code{TRUE}, the outer (redundant) \code{list} is removed.}

\item{BD.getData}{returns a list with one component per block. Each component is in turn a list
containing the stored data for the corresponding block. If \code{blocks} is of length 1 and
\code{simplify} is \code{TRUE}, the outer (redundant) \code{list} is removed.}

\item{BD.checkAndDeleteFiles}{deletes the files referenced in the input \code{bwData} if they exist.}

}

\author{
Peter Langfelder
}

\section{Warning}{The definition of \code{BlockwiseData} and the functions here 
should be considered experimental and may change in
the future.}

\seealso{
Definition of and other functions on \code{\link{BlockwiseData}}:

\code{\link{newBlockwiseData}} for creating new \code{BlockwiseData} objects;

\code{\link{mergeBlockwiseData}} for merging blockwise data structure;

\code{\link{addBlockToBlockwiseData}} for adding a new block to existing blockwise data;

}
\keyword{misc}% __ONLY ONE__ keyword per line
