\name{imputeByModule}
\alias{imputeByModule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute missing data separately in each module
}
\description{
Use \code{\link[impute]{impute.knn}} to ipmpute missing data, separately in each module. 
}
\usage{
imputeByModule(
  data, 
  labels, 
  excludeUnassigned = FALSE, 
  unassignedLabel = if (is.numeric(labels)) 0 else "grey", 
  scale = TRUE, 
  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data to be imputed, with variables (genes) in columns and observations (samples) in rows.
}
  \item{labels}{
Module labels. A vector with one entry for each column in \code{data}.
}
  \item{excludeUnassigned}{
Logical: should unassigned variables (genes) be excluded from the imputation?
}
  \item{unassignedLabel}{
The value in \code{labels} that represents unassigned variables.
}
  \item{scale}{
Logical: should \code{data} be scaled to mean 0 and variance 1 before imputation? 
}
  \item{\dots}{
Other arguments to \code{\link[impute]{impute.knn}}.
}
}
\value{
The input \code{data} with missing values imputed.
}
\author{
Peter Langfelder
}
\note{
This function is potentially faster but could give different imputed values than applying \code{impute.knn}
directly to (scaled) \code{data}.
}

\seealso{
\code{\link[impute]{impute.knn}} that does the actual imputation.
}
\keyword{misc}
