\name{numbers2colors}
\alias{numbers2colors}
\title{ Color representation for a numeric variable }
\description{
  The function creates a color represenation for the given numeric input.
}
\usage{
numbers2colors(
   x, 
   signed = NULL, 
   centered = signed, 
   lim = NULL, 
   commonLim = FALSE,
   colors = if (signed) blueWhiteRed(100) else blueWhiteRed(100)[51:100],
   naColor = "grey")
}
\arguments{
  \item{x}{ a vector or matrix of numbers. Missing values are allowed and will be assigned the color
given in \code{naColor}. If a matrix, each column of the matrix is processed separately and the return
value will be a matrix of colors. }

  \item{signed}{ logical: should \code{x} be considered signed? If \code{TRUE}, the default setting is to
use to use a palette that starts with green for the most negative values, continues with white for
values around zero and turns red for positive values. If \code{FALSE}, the default palette ranges from 
white for minimum values to red for maximum values. If not given, the behaviour is controlled by values in
\code{x}: if there are both positive and negative values, \code{signed} will be considered \code{TRUE},
otherwise \code{FALSE}.}

  \item{centered}{ logical. If \code{TRUE} and \code{signed==TRUE}, numeric value zero will
correspond to the middle of the color palette. If \code{FALSE} or \code{signed==FALSE}, the middle of
the color palette will correspond to the average of the minimum and maximum value. If neither \code{signed}
nor \code{centered} are given, \code{centered} will follow \code{signed} (see above).}

  \item{lim}{ optional specification of limits, that is numeric values that should correspond to the
first and last entry of \code{colors}. }
  
  \item{commonLim}{logical: should limits be calculated separately for each column of x, or should the
limits be the same for all columns? Only applies if \code{lim} is \code{NULL}. }

  \item{colors}{ color palette to represent the given numbers. }
  \item{naColor}{ color to represent missing values in \code{x}. }
}
\details{
  Each column of \code{x} is processed individually, meaning that the color palette is adjusted
individually for each column of \code{x}.
}
\value{
  A vector or matrix (of the same dimensions as \code{x}) of colors.
}
\author{ Peter Langfelder }

\seealso{ \code{\link{labels2colors}} for color coding of ordinal labels. }

\keyword{ misc }% __ONLY ONE__ keyword per line
