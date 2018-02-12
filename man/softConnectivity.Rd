\name{softConnectivity}
\alias{softConnectivity}
\alias{softConnectivity.fromSimilarity}
\title{ Calculates connectivity of a weighted network. }
\description{
  Given expression data or a similarity, the function constructs the adjacency matrix and for each
node calculates its connectivity, that is the sum of the adjacency to the other nodes.
}
\usage{
softConnectivity(
  datExpr, 
  corFnc = "cor", corOptions = "use = 'p'", 
  weights = NULL,
  type = "unsigned",
  power = if (type == "signed") 15 else 6, 
  blockSize = 1500, 
  minNSamples = NULL, 
  verbose = 2, indent = 0)

softConnectivity.fromSimilarity(
  similarity, 
  type = "unsigned",
  power = if (type == "signed") 15 else 6,
  blockSize = 1500, 
  verbose = 2, indent = 0)

}
\arguments{
  \item{datExpr}{ a data frame containing the expression data, with rows corresponding to samples and
columns to genes. }

  \item{similarity}{ a similarity matrix: a square symmetric matrix with entries between -1 and 1. }

  \item{corFnc}{ character string giving the correlation function to be used for the adjacency
calculation. Recommended choices are \code{"cor"} and \code{"bicor"}, but other functions can be used as
well. }

  \item{corOptions}{  character string giving further options to be passed to the correlation function. }

  \item{weights}{optional observation weights for \code{datExpr} to be used in correlation calculation.
A matrix of the same dimensions as \code{datExpr}, containing non-negative weights. Only used with Pearson
correlation.}

  \item{type}{network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. }

  \item{power}{ soft thresholding power. }

  \item{blockSize}{ block size in which adjacency is to be calculated. Too low (say below 100) may make
the calculation inefficient, while too high may cause R to run out of physical memory and slow down the
computer. Should be chosen such that an array of doubles of size (number of genes) * (block size) fits
into available physical memory.}

  \item{minNSamples}{ minimum number of samples available for the calculation of adjacency for the
adjacency to be considered valid.  If not given, defaults to the greater of \code{..minNSamples}
(currently 4) and number of samples divided by 3.  If the number of samples falls below this threshold,
the connectivity of the corresponding gene will be returned as \code{NA}. } 

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }
}
\value{
  A vector with one entry per gene giving the connectivity of each gene in the weighted network.
}
\references{ 
      Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 
}

\author{ Steve Horvath }
\seealso{ \code{\link{adjacency}} }
\keyword{ misc }
