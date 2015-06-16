                                                                     
                                                                     
                                                                     
                                             
###########################################################################
# Statistics for Microarray Analysis for R
# Discriminant analysis
#
# Date : August 21, 2000
# Last update : April 13, 2001
#
# Authors: Sandrine Dudoit, Yee Hwa (Jean) Yang, and Jane Fridlyand.
##########################################################################

##########################################################################
#                       A Red-Green Color Map
##########################################################################

########################################################################/**
#                            
# \name{rgcolors.func}
# 
# \alias{rgcolors.func}
# 
# \title{Red and Green Color Specification}
# 
# \description{
# This function creates a vector of n ``contiguous'' colors,
# corresponding to n intensities (between 0 and 1) of the red, green
# and blue primaries, with the blue intensities set to zero. The
# values returned by \code{rgcolors.func} can be used with a
# \code{col=} specification in graphics functions or in
# \code{\link{par}}.  
# }
# 
# \usage{
# rgcolors.func(n=50)
# }
# 
# \arguments{
#  \item{n}{the number of colors (>= 1) to be used in the red and
#  green palette. } 
# 
# }
# \value{a character vector of color names. Colors are specified
# directly in terms of their RGB components with a string of the form
# "\#RRGGBB", where each of the pairs RR, GG, BB consist of two
# hexadecimal digits giving a value in the range 00 to FF. 
#  }
# 
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Jane Fridlyand, \email{janef@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{plotCor}}, \code{\link{plotMat}},
# \code{\link{colors}}, \code{\link{rgb}}, \code{\link{image}}.} 
# 
# \examples{
# rgcolors.func(n=5)
# ## The following vector is returned:
# ## "#00FF00" "#40BF00" "#808000" "#BF4000" "#FF0000"
# }
# 
# \keyword{Microarray, RGB image.}
# 
#*/#######################################################################
                          
rgcolors.func<-function(n = 50) 
{
  k <- round(n/2)     
  r <- c(rep(0, k), seq(0, 1, length = k))     
  g <- c(rev(seq(0, 1, length = k)), rep(0, k))     
  res <- rgb(r, g, rep(0, 2 * k))     
  res 
}               

##########################################################################
#                Images of data matrices and correlation matrices
##########################################################################
########################################################################/**
# \name{plotCor}
# 
# \alias{plotCor}
# 
# \title{Red and Green Color Image of Correlation Matrix}
# 
# \description{
# This function produces a red and green color image of a correlation
# matrix using an RGB color specification. Increasingly positive
# correlations are represented with reds of increasing intensity, and
# increasingly negative correlations are represented with greens of
# increasing intensity.  
# }
# 
# \usage{
# plotCor(X, new=F, nrgcols=50, labels=FALSE, labcols=1, title="", ...)
# }
# 
# \arguments{
#  \item{X}{a matrix of numerical values.}
#  \item{new}{If \code{new=F}, \code{X} must already be a correlation
#  matrix. If \code{new=T}, the correlation matrix for the columns of
#  \code{X} is computed and displayed in the image.} 
#  \item{nrgcols}{the number of colors (>= 1) to be used in the red
#  and green palette.} 
#  \item{labels}{vector of character strings to be placed at the
#  tickpoints, labels for the columns of \code{X}.} 
#  \item{labcols}{colors to be used for the labels of the columns of
#  \code{X}. \code{labcols} can have either length 1, in which case
#  all the labels are displayed using the same color, or the same
#  length as \code{labels}, in which case a color is specified for the
#  label of each column of \code{X}.} 
#  \item{title}{character string, overall title for the plot.}
#  \item{\dots}{graphical parameters may also be supplied as arguments to
#           the function (see \code{\link{par}}). For comparison purposes, 
#  it is good to set \code{zlim=c(-1,1)}.}
# }
# }
# 
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{plotMat}},\code{\link{rgcolors.func}},
# \code{\link{cor.na}}, \code{\link{cor}}, \code{\link{image}},
# \code{\link{rgb}}.} 
# 
# 
# \keyword{Microarray, correlation matrix, image.}
# 
# 
#*/#######################################################################

 plotCor<-function(x, new=FALSE, nrgcols=50, labels=FALSE, labcols=1, title="", ...)
 {
#   X <- x
   n<-ncol(x)
   corr<-x
 
   if(new)
     corr<-cor(x, use = 'p')
  
   image(1:n,1:n,corr[,n:1],col=rgcolors.func(nrgcols),axes=FALSE, xlab="", ylab="",... ) 
 
  if(length(labcols)==1){
    axis(2,at=n:1,labels=labels,las=2,cex.axis=0.6,col.axis=labcols)
    axis(3,at=1:n,labels=labels,las=2,cex.axis=0.6,col.axis=labcols)
  }

  if(length(labcols)==n){
    cols<-unique(labcols)
    for(i in 1:length(cols)){
      which<-(1:n)[labcols==cols[i]]
      axis(2,at=(n:1)[which],labels=labels[which],las=2,cex.axis=0.6,col.axis=cols[i])
      axis(3,at=which,labels=labels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  mtext(title,side=3,line=3)
  box()
}

########################################################################/**
# \name{plotMat}
# 
# \alias{plotMat}
# 
# \title{Red and Green Color Image of Data Matrix}
# 
# \description{This function produces a red and green color image of a
# data matrix using an RGB color specification. Larger entries are
# represented with reds of increasing intensity, and smaller entries
# are represented with greens of increasing intensity.  
# }
# 
# \usage{
# plotMat(X, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1, ccols=1, title="",...)
# }
# 
# %- maybe also `usage' for other objects documented here.
# 
# \arguments{
#  \item{X}{a matrix of numbers.}
#  \item{nrgcols}{the number of colors (>= 1) to be used in the red
#  and green palette.} 
#  \item{rlabels}{vector of character strings to be placed at the row
#  tickpoints, labels for the rows of \code{X}.} 
#  \item{clabels}{vector of character strings to be placed at the
#  column tickpoints, labels for the columns of \code{X}.} 
#  \item{rcols}{colors to be used for the labels of the rows of
#  \code{X}. \code{rcols} can have either length 1, in which case
#  all the labels are displayed using the same color, or the same
#  length as \code{rlabels}, in which case a color is specified for the
#  label of each row of \code{X}.} 
#  \item{ccols}{colors to be used for the labels of the columns of
#  \code{X}. \code{ccols} can have either length 1, in which case
#  all the labels are displayed using the same color, or the same
#  length as \code{clabels}, in which case a color is specified for the
#  label of each column of \code{X}.} 
#  \item{title}{character string, overall title for the plot.}
#  \item{\dots}{graphical parameters may also be supplied as arguments  to
#           the function (see \code{\link{par}}).  E.g. \code{zlim=c(-3,3)}}
# }
# 
# %\references{ ~put references to the literature/web site here ~ }
# 
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{plotCor}}, \code{\link{rgcolors.func}},
# \code{\link{cor.na}}, \code{\link{cor}}, \code{\link{image}},
# \code{\link{rgb}}.} 
# 
# \examples{
# data(MouseArray)
# ##mouse.setup <- init.grid()
# ##mouse.data <- init.data() ## see \emph{init.data}
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Looking at log ratios of mouse1
# plotMat(spatial.func(mouse.lratio$M[,1], mouse.setup))
# }
# 
# \keyword{Microarray, image of data matrix.} 
# 
# 
#*/#######################################################################

plotMat<-function(x, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1, ccols=1, title="", ...)
{
#  X <-x
  n<-nrow(x)
  p<-ncol(x)	  

  image(1:p,1:n,t(x[n:1,]),col=rgcolors.func(nrgcols),axes=FALSE, xlab="", ylab="", ... ) 

  if(length(ccols)==1){
    axis(3,at=1:p,labels=clabels,las=2,cex.axis=0.6,col.axis=ccols)
      }

  if(length(ccols)==p){
    cols<-unique(ccols)
    for(i in 1:length(cols)){
      which<-(1:p)[ccols==cols[i]]
      axis(3,at=which,labels=clabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  if(length(rcols)==1){
    axis(2,at=n:1,labels=rlabels,las=2,cex.axis=0.6,col.axis=rcols)
      }

  if(length(rcols)==n){
    cols<-unique(rcols)
    for(i in 1:length(cols)){
      which<-(1:n)[rcols==cols[i]]
      axis(2,at=(n:1)[which],labels=rlabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  mtext(title,side=3,line=3)
  box()
}


