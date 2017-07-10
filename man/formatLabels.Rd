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
formatLabels(
   labels, 
   maxCharPerLine = 14, 
   maxWidth = NULL,
   maxLines = Inf,
   cex = 1,
   split = " ", 
   fixed = TRUE, 
   newsplit = split,
   keepSplitAtEOL = TRUE, 
   capitalMultiplier = 1.4,
   eol = "\n", 
   ellipsis = "...")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{labels}{Character strings to be formatted.
}
  \item{maxCharPerLine}{
Integer giving the maximum number of characters per line.
}

  \item{maxWidth}{Maximum width in user coordinates. If given, overrides \code{maxCharPerLine} above and
usually gives a much more efficient formatting.}

  \item{maxLines}{Maximum lines to retain. If a label extends past the maximum number of lines,
\code{ellipsis} is added at the end of the last line.}

  \item{cex}{Character expansion factor that the user intends to use when adding \code{labels} to the current
figure. Only used when \code{maxWidth} is specified.}

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
  \item{capitalMultiplier}{A multiplier for capital letters which typically occupy more space than lowercase
letters.}
  \item{eol}{Character string to separate lines in the output.}
  \item{ellipsis}{Chararcter string to add to the last line if the input label is longer than fits on
\code{maxLines} lines.}
}
\details{
Each given element of \code{labels} is processed independently. The character string is split using
\code{strsplit}, with \code{split} as the splitting pattern. The resulting shorter character strings are
then concatenated together with \code{newsplit} as the separator. Whenever the length (adjusted using the
capital letter multiplier) of the combined
result from the start or the previous newline character exceeds \code{maxCharPerLine}, or
\code{\link{strwidth}} exceeds \code{maxWidth}, the character specified by \code{eol}
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
