\name{pickHardThreshold}
\alias{pickHardThreshold}
\alias{pickHardThreshold.fromSimilarity}
\title{ Analysis of scale free topology for hard-thresholding. }

\description{
Analysis of scale free topology for multiple hard thresholds. The aim is to help the
user pick an appropriate threshold for network construction.
}
\usage{
pickHardThreshold(
  data, 
  dataIsExpr,
  RsquaredCut = 0.85, 
  cutVector = seq(0.1, 0.9, by = 0.05), 
  moreNetworkConcepts = FALSE,
  removeFirst = FALSE, nBreaks = 10, 
  corFnc = "cor", corOptions = "use = 'p'")

pickHardThreshold.fromSimilarity(
    similarity,
    RsquaredCut = 0.85, 
    cutVector = seq(0.1, 0.9, by = 0.05),
    moreNetworkConcepts=FALSE, 
    removeFirst = FALSE, nBreaks = 10)
}
\arguments{
  \item{data}{ expression data in a matrix or data frame. Rows correspond to samples and columns to
genes. }

  \item{dataIsExpr}{ logical: should the data be interpreted as expression (or other numeric) data, or as a
similarity matrix of network nodes? }

  \item{similarity}{ similarity matrix: a symmetric matrix with entries between -1 and 1 and unit diagonal.}

  \item{RsquaredCut}{ desired minimum scale free topology fitting index \eqn{R^2}. }

  \item{cutVector}{  a vector of hard threshold cuts for which the scale free topology fit indices
are to be calculated. }

  \item{moreNetworkConcepts}{logical: should additional network concepts be calculated? If \code{TRUE},
the function will calculate how the network density, the network heterogeneity, and the network
centralization depend on the power. For the definition of these additional network concepts, see Horvath and
Dong (2008).  PloS Comp Biol.  } 

  \item{removeFirst}{ should the first bin be removed from the connectivity histogram? }

  \item{nBreaks}{ number of bins in connectivity histograms }

  \item{corFnc}{ a character string giving the correlation function to be used in adjacency calculation.
}

  \item{corOptions}{ further options to the correlation function specified in \code{corFnc}. }

}
\details{
  The function calculates unsigned networks by thresholding the correlation matrix using thresholds given
in \code{cutVector}. For each power the scale free topology fit index is calculated
and returned along with other information on connectivity.
}
\value{
  A list with the following components:

  \item{cutEstimate}{ estimate of an appropriate hard-thresholding cut: the lowest cut for which
the scale free topology fit \eqn{R^2} exceeds \code{RsquaredCut}. If \eqn{R^2} is below \code{RsquaredCut} for all
cuts, \code{NA} is returned. }

  \item{fitIndices}{ a data frame containing the fit indices for scale free topology. The columns
contain the hard threshold, Student p-value for the correlation threshold, 
 adjusted \eqn{R^2} for the linear fit, the linear coefficient, adjusted
\eqn{R^2} for a more complicated fit models, mean connectivity, median connectivity and maximum
connectivity.  If input \code{moreNetworkConcepts} is \code{TRUE}, 3 additional columns containing network
density, centralization, and heterogeneity.}

}

\references{ Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17

Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression Network Analysis. PLoS Comput Biol
4(8): e1000117

}

\author{ Steve Horvath}

\seealso{ \code{\link{signumAdjacencyFunction}} }

\keyword{misc}
