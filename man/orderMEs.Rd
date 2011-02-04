\name{orderMEs}
\alias{orderMEs}
\title{Put close eigenvectors next to each other}
\description{
  Reorder given (eigen-)vectors such that similar ones (as measured by correlation) are next to each
other.
}
\usage{
orderMEs(MEs, greyLast = TRUE, 
         greyName = paste(moduleColor.getMEprefix(), "grey", sep=""), 
         orderBy = 1, order = NULL, 
         useSets = NULL,  verbose = 0, indent = 0)
}
\arguments{
  \item{MEs}{Module eigengenes in a multi-set format (see \code{\link{checkSets}}). A vector of lists,
with each list corresponding to one dataset and the module eigengenes in the component \code{data},
that is \code{MEs[[set]]$data[sample, module]} is the expression of the eigengene of module \code{module}
in sample
\code{sample} in dataset \code{set}. The number of samples can be different between the sets, but the
modules must be the same. }
  \item{greyLast}{Normally the color grey is reserved for unassigned genes; hence the grey module is not
a proper module and it is conventional to put it last. If this is not desired, set the parameter to
\code{FALSE}.}
  \item{greyName}{Name of the grey module eigengene.}
  \item{orderBy}{Specifies the set by which the eigengenes are to be ordered (in all other sets as well).
Defaults to the first set in \code{useSets} (or the first set, if \code{useSets} is not given).}
  \item{order}{Allows the user to specify a custom ordering.}
  \item{useSets}{Allows the user to specify for which sets the eigengene ordering is to be performed.}
  \item{verbose}{Controls verbostity of printed progress messages. 0 means silent, nonzero verbose.}
  \item{indent}{A single non-negative integer controling indentation of printed messages. 0 means no
indentation, each unit above zero adds two spaces. }
}
\details{
  Ordering module eigengenes is useful for plotting purposes. For this function the order can be
specified explicitly, or a set can be given in which the correlations of the eigengenes will determine
the order. For the latter, a hierarchical dendrogram is calculated and the order given by the dendrogram is
used for the eigengenes in all other sets.
}
\value{
  A vector of lists of the same type as \code{MEs} containing the re-ordered eigengenes.
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{\code{\link{moduleEigengenes}}, \code{\link{multiSetMEs}}, \code{\link{consensusOrderMEs}}}
\keyword{misc}
