\name{overlapTable}
\alias{overlapTable}
\title{ Calculate overlap of modules }
\description{
  The function calculates overlap counts and Fisher exact test p-values for the given two sets of module
assignments. 
}
\usage{
overlapTable(
    labels1, labels2, 
    na.rm = TRUE, ignore = NULL, 
    levels1 = NULL, levels2 = NULL)
}
\arguments{
  \item{labels1}{ a vector containing module labels. }
  \item{labels2}{ a vector containing module labels to be compared to \code{labels1}. }
  \item{na.rm}{logical: should entries missing in either \code{labels1} or \code{labels2} be removed?}
  \item{ignore}{an optional vector giving label levels that are to be ignored.}
  \item{levels1}{optional vector giving levels for \code{labels1}. Defaults to sorted unique non-missing
values in \code{labels1} that are not present in \code{ignore}.}
  \item{levels2}{optional vector giving levels for \code{labels2}. Defaults to sorted unique non-missing
values in \code{labels2} that are not present in \code{ignore}.}
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
\keyword{misc}% __ONLY ONE__ keyword per line
