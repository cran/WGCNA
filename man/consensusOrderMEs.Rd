\name{consensusOrderMEs}
\alias{consensusOrderMEs}
\title{ Put close eigenvectors next to each other in several sets. }
\description{
    Reorder given (eigen-)vectors such that similar ones (as measured by correlation) are next to each
other. This is a multi-set version of \code{\link{orderMEs}}; the dissimilarity used can be of consensus
type (for each pair of eigenvectors the consensus dissimilarity is the maximum of individual set
dissimilarities over all sets) or of majority type (for each pair of eigenvectors the consensus
dissimilarity is the average of individual set dissimilarities over all sets).
}
\usage{
consensusOrderMEs(MEs, useAbs = FALSE, useSets = NULL, 
                  greyLast = TRUE, 
                  greyName = paste(moduleColor.getMEprefix(), "grey", sep=""), 
                  method = "consensus")
}
\arguments{
  \item{MEs}{Module eigengenes of several sets in a multi-set format (see \code{\link{checkSets}}). A
vector of lists, with each list corresponding to one dataset and the module eigengenes in the component
\code{data}, that is \code{MEs[[set]]$data[sample, module]} is the expression of the eigengene of module
\code{module} in sample \code{sample} in dataset \code{set}. The number of samples can be different
between the sets, but the modules must be the same. } 
  \item{useAbs}{Controls whether vector similarity should be given by absolute value of correlation or
plain correlation.}
 \item{useSets}{Allows the user to specify for which sets the eigengene ordering is to be performed.}
  \item{greyLast}{Normally the color grey is reserved for unassigned genes; hence the grey module is not
a proper module and it is conventional to put it last. If this is not desired, set the parameter to
\code{FALSE}.}
  \item{greyName}{Name of the grey module eigengene.}
  \item{method}{A character string giving the method to be used calculating the consensus
dissimilarity. Allowed values are (abbreviations of) 
\code{"consensus"} and \code{"majority"}. The consensus dissimilarity is calculated as the
maximum of given set dissimilarities for \code{"consensus"} and as the average for \code{"majority"}.}
}
\details{
  Ordering module eigengenes is useful for plotting purposes. This function calculates the consensus
or majority
dissimilarity of given eigengenes over the sets specified by \code{useSets} (defaults to all sets).
A hierarchical dendrogram is calculated using the dissimilarity and the order given by the dendrogram is
used for the eigengenes in all other sets.
}
\value{
  A vector of lists of the same type as \code{MEs} containing the re-ordered eigengenes.
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{\code{\link{moduleEigengenes}}, \code{\link{multiSetMEs}}, \code{\link{orderMEs}}}
\keyword{misc}

