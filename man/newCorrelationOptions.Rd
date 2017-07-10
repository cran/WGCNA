\name{newCorrelationOptions}
\alias{newCorrelationOptions}
\alias{CorrelationOptions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a list of correlation options.
}
\description{
Convenience function to create a re-usable list of correlation options. 
}
\usage{
newCorrelationOptions(
      corType = c("pearson", "bicor"),
      maxPOutliers = 0.05,
      quickCor = 0,
      pearsonFallback = "individual",
      cosineCorrelation = FALSE,
      nThreads = 0,
      corFnc = if (corType=="bicor") "bicor" else "cor",
      corOptions = c(
        list(use = 'p', 
             cosine = cosineCorrelation, 
             quick = quickCor,
             nThreads = nThreads),
        if (corType=="bicor") 
           list(maxPOutliers = maxPOutliers, 
                pearsonFallback = pearsonFallback) else NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corType}{
Character specifying the type of correlation function. Currently supported options are \code{"pearson",
"bicor"}.
}
  \item{maxPOutliers}{
Maximum proportion of outliers for biweight mid-correlation. See \code{\link{bicor}}.
}
  \item{quickCor}{
Real number between 0 and 1 that controls the handling of missing data in the
calculation of correlations. See \code{\link{bicor}}.
}
  \item{pearsonFallback}{
Specifies whether the bicor calculation should revert to Pearson when median
absolute deviation (mad) is zero. Recongnized values are (abbreviations of)
\code{"none", "individual", "all"}. If set to
\code{"none"}, zero mad will result in \code{NA} for the corresponding correlation.
If set to \code{"individual"}, Pearson calculation will be used only for columns that have zero mad.
If set to \code{"all"}, the presence of a single zero mad will cause the whole variable to be treated in
Pearson correlation manner (as if the corresponding \code{robust} option was set to \code{FALSE}). 
}
  \item{cosineCorrelation}{
Logical: calculate cosine biweight midcorrelation?
Cosine bicorrelation is similar to standard bicorrelation but the median subtraction is not performed.
}
  \item{nThreads}{
A non-negative integer specifying the number of parallel threads to be used by certain
parts of correlation calculations. This option only has an effect on systems on which a POSIX thread
library is available (which currently includes Linux and Mac OSX, but excludes Windows).
If zero, the number of online processors will be used if it can be determined dynamically, otherwise
correlation calculations will use 2 threads.
}
  \item{corFnc}{
Correlation function to be called in R code. Should correspoind to the value of \code{corType} above.
}
  \item{corOptions}{
A list of options to be supplied to the correlation function (in addition to appropriate arguments \code{x}
and \code{y}).
}
}
\value{
A list containing a copy of the input arguments. The output has class \code{CorrelationOptions}.
}
\author{
Peter Langfelder
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\keyword{misc}
