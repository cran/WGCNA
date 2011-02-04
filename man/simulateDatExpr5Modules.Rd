\name{simulateDatExpr5Modules}
\alias{simulateDatExpr5Modules}
\title{ Simplified simulation of expression data}
\description{
  This function provides a simplified interface to the expression data simulation, at the cost of
considerably less flexibility.
}
\usage{
simulateDatExpr5Modules(
  nGenes = 2000,
  colorLabels = c("turquoise", "blue", "brown", "yellow", "green"),
  simulateProportions = c(0.1, 0.08, 0.06, 0.04, 0.02),
  MEturquoise, MEblue, MEbrown, MEyellow, MEgreen,
  SDnoise = 1, backgroundCor = 0.3)

}
\arguments{

  \item{nGenes}{ total number of genes to be simulated. }

  \item{colorLabels}{ labels for simulated modules. }

  \item{simulateProportions}{ a vector of length 5 giving proportions of the total number of genes to be
placed in each individual module. The entries must be positive and sum to at most 1. If the sum is less
than 1, the leftover genes will be simulated outside of modules. }

  \item{MEturquoise}{ seed module eigengene for the first module. }

  \item{MEblue}{ seed module eigengene for the second module.  }

  \item{MEbrown}{ seed module eigengene for the third module.  }

  \item{MEyellow}{ seed module eigengene for the fourth module. }

  \item{MEgreen}{ seed module eigengene for the fifth module.  }

  \item{SDnoise}{ level of noise to be added to the simulated expressions. }

  \item{backgroundCor}{ backgrond correlation. If non-zero, a component will be added to all genes such
that the average correlation of otherwise unrelated genes will be \code{backgroundCor}. }

}
\details{
  Roughly one-third of the genes are simulated with a negative correlation to their seed eigengene. See
the functions \code{\link{simulateModule}} and \code{\link{simulateDatExpr}} for more details.
}

\value{

  A list with the following components:
  \item{datExpr}{ the simulated expression data in a data frame, with rows corresponding to samples and
columns to genes. }

  \item{truemodule}{ a vector with one entry per gene containing the simulated module membership. }

  \item{datME}{a data frame containing a copy of the input module eigengenes. } 
}

\author{ Steve Horvath and Peter Langfelder }

\seealso{ 

\code{\link{simulateModule}} for simulation of individual modules;

\code{\link{simulateDatExpr}} for a more comprehensive data simulation interface. 

}
\keyword{ misc }% __ONLY ONE__ keyword per line
