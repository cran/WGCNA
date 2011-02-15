\name{labeledBarplot}
\alias{labeledBarplot}
\title{ Barplot with text or color labels. }
\description{
  Produce a barplot with extra annotation.
}
\usage{
labeledBarplot(
      Matrix, labels, 
      colorLabels = FALSE, 
      colored = TRUE, 
      setStdMargins = TRUE, 
      stdErrors = NULL, 
      cex.lab = NULL, 
      xLabelsAngle = 45,
      ...)
}
\arguments{
  \item{Matrix}{ vector or a matrix to be plotted. }

  \item{labels}{ labels to annotate the bars underneath the barplot. }

  \item{colorLabels}{ logical: should the labels be interpreted as colors? If \code{TRUE}, the bars will
be labeled by colored squares instead of text. See details.  }

  \item{colored}{ logical: should the bars be divided into segments and colored? If \code{TRUE}, assumes
the \code{labels} can be interpreted as colors, and the input
\code{Matrix} is square and the rows have the same labels as the columns. See details. }

  \item{setStdMargins}{ if \code{TRUE}, the function wil set margins \code{c(3, 3, 2, 2)+0.2}. }

  \item{stdErrors}{ if given, error bars corresponding to \code{1.96*stdErrors} will be plotted on top of
the bars. }

  \item{cex.lab}{ character expansion factor for axis labels, including the text labels underneath the
barplot. }

  \item{xLabelsAngle}{angle at which text labels under the barplot will be printed. }

  \item{\dots}{ other parameters for the function \code{\link{barplot}}. }
}
\details{

Individual bars in the barplot can be identified either by printing the text of the corresponding entry
in \code{labels} underneath the bar at the angle specified by \code{xLabelsAngle}, 
or by interpreting the \code{labels} entry as a
color (see below) and drawing a correspondingly colored square underneath the bar. 

For reasons of compatibility with other functions, \code{labels} are interpreted as colors after
stripping the first two characters from each label. For example, the label \code{"MEturquoise"} is
interpreted as the color turquoise. 

If \code{colored} is set, the code assumes that \code{labels} can be interpreted as colors, and the input
\code{Matrix} is square and the rows have the same labels as the columns. Each bar in the barplot is then
sectioned into contributions from each row entry in \code{Matrix} and is colored by the color given by the
entry in \code{labels} that corresponds to the row. 
  
}
\value{
  None.
}
\author{ Peter Langfelder }

\keyword{ hplot }% __ONLY ONE__ keyword per line
