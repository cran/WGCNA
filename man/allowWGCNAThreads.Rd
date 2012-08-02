\name{allowWGCNAThreads}
\alias{allowWGCNAThreads}
\alias{enableWGCNAThreads}
\alias{disableWGCNAThreads}
\alias{WGCNAnThreads}
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

enableWGCNAThreads(nThreads = NULL)

disableWGCNAThreads()

WGCNAnThreads()
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
\details{
\code{allowWGCNAThreads} enables parallel calculation within the compiled code in WGCNA, principally for
calculation of correlations in the presence of missing data. This function is now deprecated; use
\code{enableWGCNAThreads} instead.

\code{enableWGCNAThreads} enables parallel calculations within user-level R functions as well as within the
compiled code, and registers an
appropriate parallel calculation back-end for the operating system/platform. 

\code{disableWGCNAThreads} disables parallel processing. 

\code{WGCNAnThreads} returns the number of threads (parallel processes) that WGCNA is currently configured
to run with.
}
\value{
  \code{allowWGCNAThreads},  \code{enableWGCNAThreads}, and \code{disableWGCNAThreads} return the 
maximum number of threads WGCNA calculations will be allowed to use.
}

\note{
  Multi-threading within compiled code is not available on Windows; R code parallelization works on all
platforms. 
}

\author{
Peter Langfelder
}
\keyword{misc}
