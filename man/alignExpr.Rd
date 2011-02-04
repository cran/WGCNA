\name{alignExpr}
\alias{alignExpr}
\title{ Align expression data with given vector }
\description{
  Multiplies genes (columns) in given expression data such that their correlation with given reference
vector is non-negative.
}
\usage{
alignExpr(datExpr, y = NULL)
}
\arguments{
  \item{datExpr}{ expression data to be aligned. A data frame with columns corresponding to genes and
rows to samples. }
  \item{y}{ reference vector of length equal the number of samples (rows) in \code{datExpr} }
}
\details{
  The function basically multiplies each column in \code{datExpr} by the sign of its correlation with
\code{y}. If \code{y} is not given, the first column in \code{datExpr} will be used as the reference
vector.
}
\value{
  A data frame containing the aligned expression data, of the same dimensions as the input data frame.
}
\author{ Steve Horvath and Peter Langfelder }
\keyword{ misc }% __ONLY ONE__ keyword per line
