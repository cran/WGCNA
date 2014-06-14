\name{formatLabels}
\alias{formatLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Break long character strings into multiple lines
}
\description{
This function attempts to break lomg character strings into multiple lines by replacing a given pattern by
a newline character.
}
\usage{
formatLabels(labels, 
             maxCharPerLine = 14, 
             split = " ", fixed = TRUE, 
             newsplit = split, keepSplitAtEOL = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{labels}{Character strings to be formatted.
}
  \item{maxCharPerLine}{
Integer giving the maximum number of characters per line.
}
  \item{split}{
Pattern to be replaced by newline ('\\n') characters.
}
  \item{fixed}{
Logical: Should the pattern be interpreted literally (\code{TRUE}) or as a regular expression (\code{FALSE})?
See \code{\link{strsplit}} and its argument \code{fixed}.
}
  \item{newsplit}{
Character string to replace the occurrences of \code{split} above with.
}
  \item{keepSplitAtEOL}{
When replacing an occurrence of \code{split} with a newline character, should the \code{newsplit} be added
before the newline as well? 
}
}
\details{
Each given element of \code{labels} is processed independently. The character string is split using
\code{strsplit}, with \code{split} as the splitting pattern. The resulting shorter character strings are
then concatenated together with \code{newsplit} as the separator. Whenever the length of the combined
result from the start or the previous newline character exceeds \code{maxCharPerLine}, a newline character
is inserted (at the previous split). 

Note that individual segements (i.e., sections of the input between occurrences of \code{split}) whose
number of characters exceeds \code{maxCharPerLine} will not be split.
}
\value{
A character vector of the same length as input \code{labels}.
}
\author{
Peter Langfelder
}
\examples{
s = "A quick hare jumps over the brown fox";
formatLabels(s);
}
\keyword{misc}% __ONLY ONE__ keyword per line
