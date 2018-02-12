\name{adjacency}
\alias{adjacency}
\alias{adjacency.fromSimilarity}
\title{ Calculate network adjacency }
\description{
  Calculates (correlation or distance) network adjacency from given expression data or from a similarity.
}
\usage{
adjacency(datExpr, 
          selectCols = NULL, 
          type = "unsigned", 
          power = if (type=="distance") 1 else 6,
          corFnc = "cor", corOptions = list(use = "p"),
          weights = NULL,
          distFnc = "dist", distOptions = "method = 'euclidean'",
          weightArgNames = c("weights.x", "weights.y"))

adjacency.fromSimilarity(similarity, 
                         type = "unsigned", 
                         power = if (type=="distance") 1 else 6)

}
\arguments{
  \item{datExpr}{ data frame containing expression data. Columns correspond to genes and rows to
samples.}

  \item{similarity}{a (signed) similarity matrix: square, symmetric matrix with entries between -1 and 1. }

  \item{selectCols}{ for correlation networks only (see below); 
can be used to select genes whose adjacencies will be calculated. Should be either a
numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are
to be used. }

  \item{type}{network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}, \code{"distance"}. }

  \item{power}{soft thresholding power. }

  \item{corFnc}{ character string specifying the function to be used to calculate co-expression
similarity for correlation networks. 
Defaults to Pearson correlation. Any function returning values between -1 and 1 can be used. }

  \item{corOptions}{ character string or a list specifying additional arguments to be passed to the function given
by \code{corFnc}. Use \code{"use = 'p', method = 'spearman'"} or, equivalently, 
\code{list(use = 'p', method = 'spearman')} to obtain Spearman correlation.   }

  \item{weights}{optional observation weights for \code{datExpr} to be used in correlation calculation.
A matrix of the same dimensions as \code{datExpr}, containing non-negative weights. Only used with Pearson
correlation.}

  \item{distFnc}{ character string specifying the function to be used to calculate co-expression
similarity for distance networks. Defaults to the function \code{\link{dist}}. 
Any function returning non-negative values can be used.}

  \item{distOptions}{ character string or a list specifying additional arguments to be passed to the function given
by \code{distFnc}. For example, when the function  \code{\link{dist}} is used, the argument \code{method}
can be used to specify various ways of computing the distance. }

  \item{weightArgNames}{character vector of length 2 giving the names of the arguments to \code{corFnc} that
represent weights for variable x and y. Only used if \code{weights} are non-NULL.}

}
\details{
  The argument \code{type} determines whether a correlation (\code{type} one of \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}), or a distance network (\code{type} equal \code{"distance"}) will
be calculated. In correlation networks the adajcency is constructed from correlations (values between -1 and
1, with high numbers meaning high similarity). In distance networks, the adjacency is constructed from
distances (non-negative values, high values mean low similarity). 
  
  The function calculates the similarity of columns (genes) in \code{datExpr} by calling the function
given in \code{corFnc} (for correlation networks) or \code{distFnc} (for distance networks), 
transforms the similarity according to \code{type} and raises it to \code{power},
resulting in a weighted network adjacency matrix. If \code{selectCols} is given, the \code{corFnc} function
will be given arguments \code{(datExpr, datExpr[selectCols], ...)}; hence the returned adjacency will have
rows corresponding to all genes and columns corresponding to genes selected by \code{selectCols}.

Correlation and distance are transformed as follows: for \code{type = "unsigned"}, adjacency = |cor|^power;
for \code{type = "signed"}, adjacency = (0.5 * (1+cor) )^power; for \code{type = "signed hybrid"}, adjacency
= cor^power if cor>0 and 0 otherwise; and for \code{type = "distance"}, adjacency =
(1-(dist/max(dist))^2)^power.

The function \code{adjacency.fromSimilarity} inputs a similarity matrix, that is it skips the correlation
calculation step but is otherwise identical.

}
\value{
  Adjacency matrix of dimensions \code{ncol(datExpr)} times \code{ncol(datExpr)} (or the same dimensions 
as \code{similarity}). If \code{selectCols} was
given, the number of columns will be the length (if numeric) or sum (if boolean) of \code{selectCols}.
}

\note{When calculated from the \code{datExpr}, the network is always calculated among the columns of 
\code{datExpr} irrespective of whether a correlation or a distance network is requested. }

\references{
Bin Zhang and Steve Horvath (2005) A General Framework for Weighted Gene Co-Expression
Network Analysis, Statistical Applications in Genetics and Molecular Biology, Vol. 4 No. 1, Article 17 

Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between co-expression
modules. BMC Systems Biology 2007, 1:54
}
\author{ Peter Langfelder and Steve Horvath }
\keyword{ misc }% __ONLY ONE__ keyword per line
