\name{scaleFreeFitIndex}
\Rdversion{1.1}
\alias{scaleFreeFitIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculation of fitting statistics for evaluating scale free topology fit.
}
\description{
 The function scaleFreeFitIndex calculates several indices (fitting statistics) for evaluating scale free
topology fit.  The input is a vector (of connectivities) k. Next k is discretized into nBreaks number of
equal-width bins.  Let's denote the resulting vector dk.  The relative frequency for each bin is denoted
p.dk. 
}
\usage{
scaleFreeFitIndex(k, nBreaks = 10, removeFirst = FALSE)
}

\arguments{
  \item{k}{
  numeric vector whose components contain non-negative values
}
  \item{nBreaks}{
  positive integer. This determines the number of equal width bins.
}
  \item{removeFirst}{
  logical. If TRUE then the first bin will be removed.
}
}
\value{
  Data frame with columns
 \item{Rsquared.SFT}{the model fitting index (R.squared) from the following model lm(log.p.dk ~ log.dk)}
 \item{slope.SFT}{the slope estimate from model lm(log(p(k))~log(k))}
 \item{truncatedExponentialAdjRsquared}{the adjusted R.squared measure from the truncated exponential model
given by lm2 = lm(log.p.dk ~ log.dk + dk).}
}
\author{
Steve Horvath
}
\keyword{ misc }
