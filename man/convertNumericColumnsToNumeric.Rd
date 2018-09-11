\name{convertNumericColumnsToNumeric}
\alias{convertNumericColumnsToNumeric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert character columns that represent numbers to numeric
}
\description{
This function converts to numeric those character columns in the input that can be converted to numeric without generating
missing values except for the allowed NA representations.
}
\usage{
convertNumericColumnsToNumeric(
   data, 
   naStrings = c("NA", "NULL", "NO DATA"), 
   unFactor = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame.
}
  \item{naStrings}{
Character vector of values that are allowd to convert to \code{NA} (a missing numeric value).
}
  \item{unFactor}{
Logical: should the function first convert all factor columns to character?
}
}
\value{
A data frame with convertible columns converted to numeric.
}
\author{
Peter Langfelder
}
\keyword{misc}% __ONLY ONE__ keyword per line
