\name{plotCor}

\alias{plotCor}

\title{Red and Green Color Image of Correlation Matrix}

\description{
This function produces a red and green color image of a correlation
matrix using an RGB color specification. Increasingly positive
correlations are represented with reds of increasing intensity, and
increasingly negative correlations are represented with greens of
increasing intensity.  
}

\usage{
plotCor(x, new=FALSE, nrgcols=50, labels=FALSE, labcols=1, title="", ...)
}

\arguments{
 \item{x}{a matrix of numerical values.}
 \item{new}{If \code{new=F}, \code{x} must already be a correlation
 matrix. If \code{new=T}, the correlation matrix for the columns of
 \code{x} is computed and displayed in the image.} 
 \item{nrgcols}{the number of colors (>= 1) to be used in the red
 and green palette.} 
 \item{labels}{vector of character strings to be placed at the
 tickpoints, labels for the columns of \code{x}.} 
 \item{labcols}{colors to be used for the labels of the columns of
 \code{x}. \code{labcols} can have either length 1, in which case
 all the labels are displayed using the same color, or the same
 length as \code{labels}, in which case a color is specified for the
 label of each column of \code{x}.} 
 \item{title}{character string, overall title for the plot.}
 \item{\dots}{graphical parameters may also be supplied as arguments to
          the function (see \code{\link{par}}). For comparison purposes, 
 it is good to set \code{zlim=c(-1,1)}.}
}


\author{
  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
}

\seealso{\code{\link{plotMat}},\code{\link{rgcolors.func}},
\code{\link{cor}}, \code{\link{image}},
\code{\link{rgb}}.} 


\keyword{hplot}

