\name{newBlockwiseData}
\alias{newBlockwiseData}
\alias{BlockwiseData}
\alias{mergeBlockwiseData}
\alias{addBlockToBlockwiseData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create, merge and expand BlockwiseData objects 
}
\description{
These functions create, merge and expand BlockwiseData objects for holding in-memory or disk-backed blockwise
data. Blockwise here means
that the data is too large to be loaded or processed in one piece and is therefore split into blocks that can
be handled one by one in a divide-and-conquer manner.
}
\usage{
newBlockwiseData(
   data, 
   external = FALSE, 
   fileNames = NULL, 
   doSave = external, 
   recordAttributes = TRUE, 
   metaData = list())

mergeBlockwiseData(...)

addBlockToBlockwiseData(
   bwData,
   blockData,
   external = bwData$external,
   blockFile = NULL,
   doSave = external,
   recordAttributes = !is.null(bwData$attributes),
   metaData = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A list in which each component carries the data of a single block.
}
  \item{external}{
Logical: should the data be disk-backed (\code{TRUE}) or in-memory (\code{FALSE})?
}
  \item{fileNames}{
When \code{external} is \code{TRUE}, this argument must be a 
character vector of the same length as \code{data}, giving the file names for the data to be saved to, or
where the data is already located.
}
  \item{doSave}{
Logical: should data be saved? If this is \code{FALSE}, it is the user's responsibility to ensure the files
supplied in \code{fileNames} already exist and contain the expected data.
}
  \item{recordAttributes}{
Logical: should \code{attributes} of the given data be recorded within the object?
}
  \item{metaData}{
A list giving any additional meta-data for \code{data} that should be attached to the object.
}
  \item{bwData}{An existing \code{BlockwiseData} object.}

  \item{blockData}{A vector, matrix or array carrying the data of a single block.}

  \item{blockFile}{
File name where data contained in \code{blockData} should be saved.
}
  \item{...}{One or more objects of class \code{BlockwiseData}.}
}
\details{
Several functions in this package use the concept of blockwise, or "divide-and-conquer", analysis. The
BlockwiseData class is meant to hold the blockwise data, or all necessary information about blockwise data
that is saved in disk files. 

The data can be stored in disk files (one file per block) or in-memory. In memory storage is provided so that
same code can be used for both smaller (single-block) data where disk storage could slow down operations as
well as larger data sets where disk storage and block by block analysis are necessary.
}
\value{
All three functions return a list with the class set to \code{"BlockwiseData"}, containing the following components:
\item{external}{Copy of the input argument \code{external}}
\item{data}{If \code{external} is \code{TRUE}, an empty list, otherwise a copy of the input \code{data}.}
\item{fileNames}{Copy of the input argument \code{fileNames}.}
\item{lengths}{A vector of lengths (results of \code{\link{length}}) of elements of \code{data}.}
\item{attributes}{If input \code{recordAttributes} is \code{TRUE}, a list with one component per block
(component of \code{data}); each component is in turn a list of attributes of that component of \code{data}.}
\item{metaData}{A copy of the input \code{metaData}.}
}
\author{
Peter Langfelder
}

\section{Warning}{The definition of \code{BlockwiseData} should be considered experimental and may change in
the future.}

\seealso{
Other functions on \code{BlockwiseData}: 

\code{\link{BD.getData}} for retrieving data

\code{\link{BD.actualFileNames}} for retrieving file names of files containing data;

\code{\link{BD.nBlocks}} for retrieving the number of blocks;

\code{\link{BD.blockLengths}} for retrieving block lengths;

\code{\link{BD.getMetaData}} for retrieving metadata;

\code{\link{BD.checkAndDeleteFiles}} for deleting files of an unneeded object.
}
\keyword{misc}
