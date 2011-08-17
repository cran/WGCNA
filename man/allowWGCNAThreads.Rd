\name{allowWGCNAThreads}
\alias{allowWGCNAThreads}
\alias{disableWGCNAThreads}
\title{
Allow and disable multi-threading for certain WGCNA calculations
}
\description{
These functions allow and disable multi-threading for WGCNA calculations that can optionally be
multi-threaded, which includes all functions using \code{\link{cor}}
or \code{\link{bicor}} functions. 
}
\usage{
allowWGCNAThreads(nThreads = NULL)

disableWGCNAThreads()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nThreads}{
Number of threads to allow. If not given, the number of processors online (as reported by system
configuration) will be used. There appear to be some cases where the automatically-determined number is
wrong; please check the output to see that the number of threads makes sense. Except for testing and/or
torturing your system, the number of threads should be no more than the number of actual processors/cores.
}
}
\value{
  Maximum number of threads WGCNA calculations will be allowed to use.
}

\note{
  Multi-threading within WGCNA is not available on Windows.
}

\author{
Peter Langfelder
}
\keyword{misc}
