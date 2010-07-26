\name{fixDataStructure}
\alias{fixDataStructure}
\title{Put single-set data into a form useful for multiset calculations. }
\description{
  Encapsulates single-set data in a wrapper that makes the data suitable for functions working on
multiset data collections.
}
\usage{
fixDataStructure(data, verbose = 0, indent = 0)
}
\arguments{
  \item{data}{ A dataframe, matrix or array with two dimensions to be encapsulated. }
  \item{verbose}{Controls verbosity. 0 is silent. }
  \item{indent}{Controls indentation of printed progress messages. 0 means no indentation, every unit
adds two spaces.}
}
\details{
  For multiset calculations, many quantities (such as expression data, traits, module eigengenes etc) are
presented by a common structure, a vector of lists (one list for each set) where each list has a
component \code{data} that contains the actual (expression, trait, eigengene) data for the corresponding
set in the form of a dataframe. This funtion creates a vector of lists of length 1 and fills the
component \code{data} with the content of parameter \code{data}.
}
\value{
  As described above, input data in a format suitable for functions operating on multiset data
collections.
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{ \code{\link{checkSets}}}
\examples{

singleSetData = matrix(rnorm(100), 10,10);
encapsData = fixDataStructure(singleSetData);
length(encapsData)
names(encapsData[[1]])
dim(encapsData[[1]]$data)
all.equal(encapsData[[1]]$data, singleSetData);

}

\keyword{misc}
