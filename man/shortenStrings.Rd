\name{shortenStrings}
\alias{shortenStrings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shorten given character strings by truncating at a suitable separator.
}
\description{
This function shortens given character strings so they are not longer than a given maximum length.
}
\usage{
shortenStrings(strings, maxLength = 25, minLength = 10, 
              split = " ", fixed = TRUE,
              ellipsis = "...", countEllipsisInLength = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strings}{
Character strings to be shortened.
}
  \item{maxLength}{
Maximum length (number of characters) in the strings to be retained. See details for when the returned
strings can exceed this length.
}
  \item{minLength}{
Minimum length of the returned strings. See details.
}
  \item{split}{
Character string giving the split at which the strings can be truncated. This can be a literal string or a 
regular expression (if the latter, \code{fixed} below must be set to \code{FALSE}).
}
  \item{fixed}{
Logical: should \code{split} be interpreted as a literal specification (\code{TRUE}) or as a regular
expression (\code{FALSE})? }

  \item{ellipsis}{
   Character string that will be appended to every shorten string, to indicate that the string has been
shortened. }

  \item{countEllipsisInLength}{
   Logical: should the length of the ellipsis count toward the minimum and maximum length?}
}
\details{
Strings whose length (number of characters) is at most \code{maxLength} are returned unchanged. For those
that are longer, the function uses \code{\link{gregexpr}} to search for the occurrences of \code{split} in
each given 
character string. If such occurrences are found at positions between \code{minLength} and \code{maxLength},
the string will be truncated at the last such \code{split}; otherwise, the string will be truncated at
\code{maxLength}. The \code{ellipsis} is appended to each truncated string.
}
\value{
A character vector of strings, shortened as necessary. If the input \code{strings} had non-NULL dimensions
and dimnames, these are copied to the output.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{gregexpr}}, the workhorse pattern matching function
\code{\link{formatLabels}} for splitting strings into multiple lines
}
\keyword{ misc }% __ONLY ONE__ keyword per line
