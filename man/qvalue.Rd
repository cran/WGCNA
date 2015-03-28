\name{qvalue}
\alias{qvalue}
\title{Estimate the q-values for a given set of p-values}
\description{
  Estimate the q-values for a given set of p-values.  The q-value of a
  test measures the proportion of false positives incurred (called the
  false discovery rate) when that particular test is called significant.
}
\usage{
qvalue(p, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=NULL, robust=FALSE,
  smooth.df=3, smooth.log.pi0=FALSE)
}
\arguments{
  \item{p}{A vector of p-values (only necessary input)}
  \item{lambda}{The value of the tuning parameter to estimate
    \eqn{\pi_0}{pi_0}. Must be in [0,1). Optional, see Storey (2002).}
  \item{pi0.method}{Either "smoother" or "bootstrap"; the method for
    automatically choosing tuning parameter in the estimation of \eqn{\pi_0}{pi_0}, 
    the proportion of true null hypotheses}
  \item{fdr.level}{A level at which to control the FDR. Must be in (0,1]. Optional; if this is
    selected, a vector of TRUE and FALSE is returned that specifies
    whether each q-value is less than fdr.level or not.}
  \item{robust}{An indicator of whether it is desired to make the
    estimate more robust for small p-values and a direct finite sample 
    estimate of pFDR. Optional.}
  \item{smooth.df}{Number of degrees-of-freedom to use when estimating \eqn{\pi_0}{pi_0} 
    with a smoother. Optional.}
  \item{smooth.log.pi0}{If TRUE and \code{pi0.method} = "smoother", \eqn{\pi_0}{pi_0} will be 
    estimated by applying a smoother to a scatterplot of \eqn{log} \eqn{\pi_0}{pi_0} estimates 
    against the tuning parameter \eqn{\lambda}{lambda}. Optional.}
}
\details{
  If no options are selected, then the method used to estimate \eqn{\pi_0}{pi_0} is
  the smoother method described in Storey and Tibshirani (2003). The
  bootstrap method is described in Storey, Taylor & Siegmund (2004).
}
\value{
  A list containing:
\item{call}{function call}
\item{pi0}{an estimate of the proportion of null p-values}
\item{qvalues}{a vector of the estimated q-values (the main quantity of
  interest)} 
\item{pvalues}{a vector of the original p-values}
\item{significant}{if fdr.level is specified, and indicator of whether the
  q-value fell below fdr.level (taking all such q-values to be significant
  controls FDR at level fdr.level)}
}
\note{This function is adapted from package qvalue. The reason we provide our own copy is that 
package qvalue contains additional functionality that relies on Tcl/Tk which has led to multiple problems.
Our copy does not require Tcl/Tk.}

\references{
  Storey JD. (2002) A direct approach to false discovery rates. Journal
  of the Royal Statistical Society, Series B, 64: 479-498.
  
  Storey JD and Tibshirani R. (2003) Statistical significance for
  genome-wide experiments. Proceedings of the National Academy of Sciences, 
  100: 9440-9445. 

  Storey JD. (2003) The positive false discovery rate: A Bayesian
  interpretation and the q-value. Annals of Statistics, 31: 2013-2035.  

  Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
  conservative point estimation, and simultaneous conservative
  consistency of false discovery rates: A unified approach. Journal of
  the Royal Statistical Society, Series B, 66: 187-205.
}
\author{John D. Storey \email{jstorey@u.washington.edu}, adapted for WGCNA by Peter Langfelder}
\keyword{misc}
