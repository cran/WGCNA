\name{pickSoftThreshold}
\alias{pickSoftThreshold}
\alias{pickSoftThreshold.fromSimilarity}
\title{ Analysis of scale free topology for soft-thresholding }
\description{
  Analysis of scale free topology for multiple soft thresholding powers. The aim is to help the user pick
an appropriate soft-thresholding power for network construction.
}
\usage{
pickSoftThreshold(
  data, 
  dataIsExpr = TRUE,
  weights = NULL,
  RsquaredCut = 0.85, 
  powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
  removeFirst = FALSE, nBreaks = 10, blockSize = NULL, 
  corFnc = cor, corOptions = list(use = 'p'), 
  networkType = "unsigned",
  moreNetworkConcepts = FALSE,
  gcInterval = NULL,
  verbose = 0, indent = 0)

pickSoftThreshold.fromSimilarity(
    similarity,
    RsquaredCut = 0.85, 
    powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
    removeFirst = FALSE, nBreaks = 10, blockSize = 1000,
    networkType = "unsigned", 
    moreNetworkConcepts=FALSE, 
    verbose = 0, indent = 0)

}
\arguments{
  \item{data}{ expression data in a matrix or data frame. Rows correspond to samples and columns to
genes. }

  \item{dataIsExpr}{ logical: should the data be interpreted as expression (or other numeric) data, or as a
similarity matrix of network nodes? }

    \item{weights}{optional observation weights for \code{data} to be used in correlation calculation.
A matrix of the same dimensions as \code{datExpr}, containing non-negative weights. Only used with Pearson
correlation.}

  \item{similarity}{ similarity matrix: a symmetric matrix with entries between -1 and 1 and unit diagonal.}

  \item{RsquaredCut}{ desired minimum scale free topology fitting index \eqn{R^2}. }

  \item{powerVector}{ a vector of soft thresholding powers for which the scale free topology fit indices
are to be calculated. }

  \item{removeFirst}{ should the first bin be removed from the connectivity histogram? }

  \item{nBreaks}{ number of bins in connectivity histograms }

  \item{blockSize}{ block size into which the calculation of connectivity should be broken up. If not given,
a suitable value will be calculated using function \code{blockSize} and printed if \code{verbose>0}. If R runs
into memory problems, decrease this value. }

  \item{corFnc}{the correlation function to be used in adjacency calculation. }

  \item{corOptions}{ a list giving further options to the correlation function specified in \code{corFnc}. }

  \item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link{adjacency}}. }

  \item{moreNetworkConcepts}{logical: should additional network concepts be calculated? If \code{TRUE}, the
function will calculate how the network density, the network heterogeneity, and the network centralization
depend on the power. For the definition of these additional network concepts, see Horvath and Dong (2008).
PloS Comp Biol.  }

  \item{gcInterval}{a number specifying in interval (in terms of individual genes) in which garbage
collection will be performed. The actual interval will never be less than \code{blockSize}.}

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }
  
  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\details{
  The function calculates weighted networks either by
interpreting \code{data} directly as similarity, or first transforming it to similarity of the type
specified by \code{networkType}.
The weighted networks are obtained by raising the similarity to the powers given in \code{powerVector}. 
For each power the scale free topology fit index is calculated
and returned along with other information on connectivity. 

On systems with multiple cores or processors, 
the function pickSoftThreshold takes advantage of parallel processing if the function 
\code{\link{enableWGCNAThreads}} has been called to allow parallel processing and set up the parallel
calculation back-end.
}
\value{
  A list with the following components:

  \item{powerEstimate}{ estimate of an appropriate soft-thresholding power: the lowest power for which
the scale free topology fit \eqn{R^2} exceeds \code{RsquaredCut}. If \eqn{R^2} is below \code{RsquaredCut} 
for all powers, \code{NA} is returned. }

  \item{fitIndices}{ a data frame containing the fit indices for scale free topology. The columns
contain the soft-thresholding power, adjusted \eqn{R^2} for the linear fit, the linear coefficient, adjusted
\eqn{R^2} for a more complicated fit models, mean connectivity, median connectivity and maximum
connectivity. If input \code{moreNetworkConcepts} is \code{TRUE}, 3 additional columns containing network
density, centralization, and heterogeneity.} 

}
\references{ Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17

Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression Network Analysis. PLoS Comput Biol
4(8): e1000117 
}

\author{ Steve Horvath and Peter Langfelder }

\seealso{ \code{\link{adjacency}}, \code{\link{softConnectivity}} }

\keyword{misc}

