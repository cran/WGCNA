\name{factorizeNonNumericColumns}
\alias{factorizeNonNumericColumns}
\title{
Turn non-numeric columns into factors
}
\description{
Given a data frame, this function turns non-numeric columns into factors.
}
\usage{
factorizeNonNumericColumns(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame. Non-data frame inputs (e.g., a matrix) are coerced to a data frame.
}
}
\details{
A column is considered numeric if its storage mode is numeric or if it is a character vector, it only contains
character representations of numbers and possibly missing values encoded as "NA", "NULL", "NO DATA".
}
\value{
The input data frame with non-numeric columns turned into factors.
}
\author{
Peter Langfelder
}
\keyword{misc}% __ONLY ONE__ keyword per line
