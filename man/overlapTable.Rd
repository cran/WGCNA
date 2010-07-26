\name{overlapTable}
\alias{overlapTable}
\title{ Calculate overlap of modules }
\description{
  The function calculates overlap counts and Fisher exact test p-values for the given two sets of module
assignments. 
}
\usage{
overlapTable(labels1, labels2)
}
\arguments{
  \item{labels1}{ a vector containing module labels. }
  \item{labels2}{ a vector containing module labels to be compared to \code{labels1}. }
}
\value{
  A list with the following components:
  \item{countTable}{a matrix whose rows correspond to modules (unique labels) in \code{labels1} and whose
columns correspond to modules (unique labels) in \code{labels2}, giving the number of objects in the
intersection of the two respective modules. }

  \item{pTable}{a matrix whose rows correspond to modules (unique labels) in \code{labels1} and whose
columns correspond to modules (unique labels) in \code{labels2}, giving Fisher's exact test 
significance p-values for the overlap of the two respective modules. }

}
\author{ Peter Langfelder }
\seealso{ \code{\link{fisher.test}}, \code{\link{matchLabels}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
