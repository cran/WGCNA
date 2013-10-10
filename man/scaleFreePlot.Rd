\name{scaleFreePlot}
\alias{scaleFreePlot}
\title{ Visual check of scale-free topology }
\description{
  A simple visula check of scale-free network ropology. 
}
\usage{
scaleFreePlot(
  connectivity, 
  nBreaks = 10, 
  truncated = FALSE,  
  removeFirst = FALSE, 
  main = "", ...)
}
\arguments{
  \item{connectivity}{ vector containing network connectivities. }

  \item{nBreaks}{ number of breaks in the connectivity dendrogram. }

  \item{truncated}{ logical: should a truncated exponential fit be calculated and plotted in addition to
the linear one? }

  \item{removeFirst}{ logical: should the first bin be removed from the fit? }

  \item{main}{ main title for the plot. }

  \item{\dots}{ other graphical parameter to the \code{plot} function. }
}
\details{
  The function plots a log-log plot of a histogram of the given \code{connectivities}, and fits a linear
model plus optionally a truncated exponential model. The \eqn{R^2} of the fit can be considered an index
of the scale freedom of the network topology. 
}
\value{
  None.
}
\references{ 
      Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 
}

\author{ Steve Horvath }
\seealso{ \code{\link{softConnectivity}} for connectivity calculation in weigheted networks. }
\keyword{ misc }% __ONLY ONE__ keyword per line
