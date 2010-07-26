\name{correlationPreservation}
\alias{correlationPreservation}
\title{ Preservation of eigengene correlations }
\description{
  Calculates a summary measure of  preservation of eigengene correlations across data sets
}
\usage{
correlationPreservation(multiME, setLabels, excludeGrey = TRUE, greyLabel = "grey")
}
\arguments{
  \item{multiME}{consensus module eigengenes in a multi-set format. A vector of lists with one list
corresponding to each set. Each list must contain a component \code{data} that is a data frame whose
columns are consensus module eigengenes. } 
  \item{setLabels}{names to be used for the sets represented in \code{multiME}.}
  \item{excludeGrey}{logical: exclude the 'grey' eigengene from preservation measure?}
  \item{greyLabel}{module label corresponding to the 'grey' module. Usually this will be the
character string \code{"grey"} if the labels are colors, and the number 0 if the labels are numeric.}
}
\details{
  The function calculates the preservation of correlation of each eigengene with all other
eigengenes (optionally except the 'grey' eigengene) in all pairs of sets.  
}
\value{
  A data frame whose rows correspond to consensus module eigengenes given in the input
\code{multiME}, and columns correspond to all possible set comparisons. The two sets compared in
each column are indicated in the column name.
}
\references{ Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships
between co-expression modules. BMC Systems Biology 2007, 1:54}
\author{ Peter Langfelder }
\seealso{ \code{\link{multiSetMEs}} and module\code{\link{checkSets}} in package moduleColor for
more on eigengenes and the multi-set format }
\keyword{ misc }
