\name{checkSets}
\alias{checkSets}
\title{Check structure and retrieve sizes of a group of datasets.  }
\description{
  Checks whether given sets have the correct format and retrieves dimensions.
}
\usage{
checkSets(data, checkStructure = FALSE, useSets = NULL)
}
\arguments{
  \item{data}{ A vector of lists; in each list there must be a component named \code{data} whose content
is a matrix or dataframe or array of dimension 2. }
  \item{checkStructure}{If \code{FALSE}, incorrect structure of \code{data} will trigger an error. If
\code{TRUE}, an appropriate flag (see output) will be set to indicate whether \code{data} has correct
structure.}
  \item{useSets}{Optional specification of entries of the vector \code{data} that are to be checked.
Defaults to all components. This may be useful when \code{data} only contains information for some
of the sets.}
}
\details{
  For multiset calculations, many quantities (such as expression data, traits, module eigengenes etc) are
presented by a common structure, a vector of lists (one list for each set) where each list has a
component \code{data} that contains the actual (expression, trait, eigengene) data for the corresponding
set in the form of a dataframe. This funtion checks whether \code{data} conforms to this convention and
retrieves some basic dimension information (see output).
}
\value{ A list with components
  \item{nSets}{Number of sets (length of the vector \code{data}).}
  \item{nGenes}{Number of columns in the \code{data} components in the lists. This number must be the
same for all sets.}
  \item{nSamples}{A vector of length \code{nSets} giving the number of rows in the \code{data}
components.}
  \item{structureOK}{Only set if the argument \code{checkStructure} equals \code{TRUE}. 
            The value is \code{TRUE} if the paramter \code{data} passes a few tests of its
structure, and \code{FALSE} otherwise. The tests are not exhaustive and are meant to catch obvious user
errors rather than be bulletproof.}
}

\author{ 
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
