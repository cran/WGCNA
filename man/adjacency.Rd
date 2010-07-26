\name{adjacency}
\alias{adjacency}
\title{ Calculate network adjacency }
\description{
  Calculates network adjacency from given expression data.
}
\usage{
adjacency(datExpr, selectCols = NULL, power = 6, type = "unsigned", corFnc = "cor", corOptions = "use = 'p'")
}
\arguments{
  \item{datExpr}{ data frame containing expression data. Columns correspond to genes and rows to
samples.}

  \item{selectCols}{ can be used to select genes whose adjacencies will be calculated. Should be either a
numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are
to be used. }

  \item{power}{soft thresholding power. }

  \item{type}{network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. }

  \item{corFnc}{ character string specifying the function to be used to calculate co-expression
similarity. Defaults to Pearson correlation. Any function returning values between -1 and 1 can be used. }

  \item{corOptions}{ character string specifying additional arguments to be passed to the function given
by \code{corFnc}. Use \code{"use = 'p', method = 'Spearman'"} to obtain Spearman correlation.   }

}
\details{
  The function calculates the similarity of columns (genes) in \code{datExpr} by calling the function
given in corFnc, transforms the similarity according to \code{type} and raises it to \code{power},
resulting in a weighted network adjacency matrix. If \code{selectCols} is given, the \code{corFnc} function
will be given arguments \code{(datExpr, datExpr[selectCols], ...)}; hence the returned adjacency will have
rows corresponding to all genes and columns corresponding to genes selected by \code{selectCols}.
}
\value{
  Adjacency matrix of dimensions \code{nrow(datExpr)} times \code{nrow(datExpr)}. If \code{selectCols} was
given, the number of columns will be the length (if numeric) or sum (if boolean) of \code{selectCols}.
}
\references{
Bin Zhang and Steve Horvath (2005) A General Framework for Weighted Gene Co-Expression
Network Analysis, Statistical Applications in Genetics and Molecular Biology, Vol. 4 No. 1, Article 17 

Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between co-expression
modules. BMC Systems Biology 2007, 1:54
}
\author{ Peter Langfelder and Steve Horvath }
\keyword{ misc }% __ONLY ONE__ keyword per line
