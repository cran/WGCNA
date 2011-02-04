###########################################################################
# Statistics for Microarray Analysis
# Function dealing with NA's
#
# Date : March 19, 2001
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
#
# Feb 2, 2004 - add "..." to sum.na and prod.na
#
##########################################################################

##########################################################################
# Basic statistics functions that are able to handle missing values
##########################################################################

########################################################################/**
# \name{na}
# 
# \alias{log.na}
# \alias{sum.na}
# \alias{mean.na}
# \alias{var.na}
# \alias{cor.na}
# \alias{quantile.na}
# \alias{length.na}
# \alias{order.na}
# \alias{scale.na}
# \alias{prod.na}
# 
# \title{Basic Statistical Functions for Handling Missing Values}
# 
# \description{
# Basic statistical functions for handling missing values or NA. \cr 
# In \code{log.na}, \code{sum.na}, \code{mean.na} and \code{var.na},
# \code{quantile.na}, \code{length.na}, missing values are omitted
# from the calculation. \cr 
# The function \code{cor.na} calls \code{cor} with the argument
# \code{use="pairwise.complete.obs"}. \cr 
# The function \code{order.na} only handles vector arguments and not
# lists.  However, it gives the option of omitting the NAs
# (\code{na.last=NA}), of placing the NAs at the start of the ordered
# vector (\code{na.last=F}) or at the end (\code{na.last=T}). \cr 
# The function \code{scale.na} is a modified version of
# \code{\link{scale}} which allows NAs in the variance calculation. If
# \code{scale = T}, the function \code{f} in \code{scale.na} uses
# \code{var.na} to perform the variance calculation.
# The function \code{prod.na} is similar to the \code{\link{prod}}
# function with \code{na.rm=TRUE}. This function returns the product of
# all the values present in its arguments, omitting any missing values.
# }
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{log}}, \code{\link{sum}}, \code{\link{mean}},
#   \code{\link{var}}, \code{\link{cor}}, \code{\link{order}},
#   \code{\link{scale}}, \code{link{prod}}.}
# 
# \keyword{log, sum, mean, variance, correlation, order, scale,
# product, missing values, NA.} 
# 
#*/#########################################################################

 
mean.na <- function(x,...)
{
        mean(x[!(is.na(x) | is.infinite(x))])
}
 
 var.na <- function(x)
{
        res <- NA
        tmp <- !(is.na(x) | is.infinite(x))
        if(sum(tmp) > 1)
                res <- var(x[tmp])
        res
}

cor.na <- function(x)
{
  cor(x, use="pairwise.complete.obs")
}

sum.na <- function(x,...)
{
        res <- NA
        tmp <- !(is.na(x) | is.infinite(x))
        if(sum(tmp) > 0)
                res <- sum(x[tmp])
        res
}


length.na <- function(x, ...)
{
   tmp <- !(is.na(x) | is.infinite(x))
   length(x[tmp],...)
 }

log.na <- function(x, ...)
{
  log(ifelse(x > 0, x, NA), ...)
}


quantile.na <- function(x, ...)
 {          
   tmp <- !(is.na(x) | is.infinite(x))
   quantile(x[tmp],...)
 }
   
order.na <- function (x, na.last = TRUE) 
{
    y <- order(x)
    n <- sum(is.na(x))
    tmp <- (length(x) - n + 1):length(x)
    if (!is.na(na.last)) {
        if (na.last) 
            res <- y
        if (!na.last)
          {
            if(n == 0)
              res <- y
            else
              res <- c(y[tmp], y[-tmp])
          }
      }
    if (is.na(na.last)) {
        warning("NA's discarded")
        res <- y[-tmp]
    }
    res
}

scale.na<-function(x, center = TRUE, scale = TRUE)
{
  x <- as.matrix(x)
  nc <- ncol(x)

  if (is.logical(center)) {
     if (center)
       x <- sweep(x, 2, apply(x, 2, mean, na.rm=TRUE))
    }
  else if (is.numeric(center) && (length(center) == nc))
    x <- sweep(x, 2, center)
  else
    stop("Length of center must equal the number of columns of x")
  
  if (is.logical(scale)) {
    if (scale) {
      f <- function(v) {
        sqrt(var.na(v))
      }
      x <- sweep(x, 2, apply(x, 2, f), "/")                   
    }
    }
  else if (is.numeric(scale) && length(scale) == nc)
    x <- sweep(x, 2, scale, "/")
  else
    stop("Length of scale must equal the number of columns of x")
    x
}

prod.na <- function (x,...) 
{
  prod(x[!(is.na(x) | is.infinite(x))])
}


##########################################################################
#                                End of file
##########################################################################
