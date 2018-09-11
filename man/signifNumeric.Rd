\name{signifNumeric}
\alias{signifNumeric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Round numeric columns to given significant digits.
}
\description{
This function applies \code{link{signif}} (or possibly other rounding function) to numeric, non-integer
columns of a given data frame.
}
\usage{
signifNumeric(x, digits, fnc = "signif")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Input data frame, matrix or matrix-like object that can be coerced to a data frame.
}
  \item{digits}{
Significant digits to retain.
}
  \item{fnc}{
The rounding function. Typically either \code{\link{signif}} or \code{\link{round}}.
}
}
\details{
The function \code{fnc} is applied to each numeric column that contains at least one non-integer (i.e., at
least one element that does not equal its own \code{round}).
}
\value{
The transformed data frame.
}
\author{
Peter Langfelder
}
\seealso{
The rounding functions \code{\link{signif}} and \code{\link{round}}.
}
\examples{
  df = data.frame(text = letters[1:3], ints = c(1:3)+234, nonints = c(0:2) + 0.02345);
  df;
  signifNumeric(df, 2);
  signifNumeric(df, 2, fnc = "round");
}
\keyword{misc}% __ONLY ONE__ keyword per line
