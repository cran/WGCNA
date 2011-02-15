\name{propVarExplained}
\alias{propVarExplained}
\title{ Proportion of variance explained by eigengenes. }
\description{
  This function calculates the proportion of variance of genes in each module explained by the respective
module eigengene. 
}
\usage{
propVarExplained(datExpr, colors, MEs, corFnc = "cor", corOptions = "use = 'p'")
}
\arguments{

 \item{datExpr}{  expression data. A data frame in which columns are genes and rows ar samples. NAs are
allowed and will be ignored. }

  \item{colors}{ a vector giving module assignment for genes given in \code{datExpr}. Unique values
should correspond to the names of the eigengenes in \code{MEs}. }

  \item{MEs}{ a data frame of module eigengenes in which each column is an eigengene and each row
corresponds to a sample.  }

  \item{corFnc}{ character string containing the name of the function to calculate correlation. Suggested
functions include \code{"cor"} and \code{"bicor"}. }

  \item{corOptions}{ further argument to the correlation function. }

}
\details{
  For compatibility with other functions, entries in \code{color} are matched to a substring of
\code{names(MEs)} starting at position 3. For example, the entry \code{"turquoise"} in \code{colors} will
be matched to the eigengene named \code{"MEturquoise"}. The first two characters of the eigengene name
are ignored and can be arbitrary.
}
\value{
  A vector with one entry per eigengene containing the proportion of variance of the module explained by
the eigengene. 
}
\author{ Peter Langfelder }
\seealso{ \code{\link{moduleEigengenes}} }
\keyword{ misc }
