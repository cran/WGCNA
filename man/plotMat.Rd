\name{plotMat}

\alias{plotMat}

\title{Red and Green Color Image of Data Matrix}

\description{This function produces a red and green color image of a
data matrix using an RGB color specification. Larger entries are
represented with reds of increasing intensity, and smaller entries
are represented with greens of increasing intensity.  
}

\usage{
plotMat(x, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1, ccols=1, title="",...)
}

%- maybe also `usage' for other objects documented here.

\arguments{
 \item{x}{a matrix of numbers.}
 \item{nrgcols}{the number of colors (>= 1) to be used in the red
 and green palette.} 
 \item{rlabels}{vector of character strings to be placed at the row
 tickpoints, labels for the rows of \code{x}.} 
 \item{clabels}{vector of character strings to be placed at the
 column tickpoints, labels for the columns of \code{x}.} 
 \item{rcols}{colors to be used for the labels of the rows of
 \code{x}. \code{rcols} can have either length 1, in which case
 all the labels are displayed using the same color, or the same
 length as \code{rlabels}, in which case a color is specified for the
 label of each row of \code{x}.} 
 \item{ccols}{colors to be used for the labels of the columns of
 \code{x}. \code{ccols} can have either length 1, in which case
 all the labels are displayed using the same color, or the same
 length as \code{clabels}, in which case a color is specified for the
 label of each column of \code{x}.} 
 \item{title}{character string, overall title for the plot.}
 \item{\dots}{graphical parameters may also be supplied as arguments  to
          the function (see \code{\link{par}}).  E.g. \code{zlim=c(-3,3)}}
}

%\references{ ~put references to the literature/web site here ~ }


\author{
  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
}

\seealso{\code{\link{plotCor}}, \code{\link{rgcolors.func}},
\code{\link{cor}}, \code{\link{image}},
\code{\link{rgb}}.} 

\keyword{hplot}


