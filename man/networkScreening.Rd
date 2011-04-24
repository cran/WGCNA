\name{networkScreening}
\alias{networkScreening}
\title{ Identification of genes related to a trait }
\description{
  This function blends standard and network approaches to selecting genes (or variables in general) highly
related to a given trait.
}
\usage{
networkScreening(y, datME, datExpr, 
                 corFnc = "cor", corOptions = "use = 'p'",
                 oddPower = 3, 
                 blockSize = 1000, 
                 minimumSampleSize = ..minNSamples,
                 addMEy = TRUE, removeDiag = FALSE, 
                 weightESy = 0.5, getQValues = TRUE)
}
\arguments{
  \item{y}{ clinical trait given as a numeric vector (one value per sample) }

  \item{datME}{ data frame of module eigengenes }

  \item{datExpr}{ data frame of expression data }

  \item{corFnc}{ character string specifying the function to be used to calculate co-expression
similarity. Defaults to Pearson correlation. Any function returning values between -1 and 1 can be used. }

  \item{corOptions}{ character string specifying additional arguments to be passed to the function given
by \code{corFnc}. Use \code{"use = 'p', method = 'spearman'"} to obtain Spearman correlation.   }

  \item{oddPower}{ odd integer used as a power to raise module memberships and significances }

  \item{blockSize}{ block size to use for calculations with large data sets }

  \item{minimumSampleSize}{ minimum acceptable number of samples. Defaults to the default minimum number of
samples used throughout the WGCNA package, currently 4.}

  \item{addMEy}{ logical: should the trait be used as an additional "module eigengene"?}

  \item{removeDiag}{ logical: remove the diagonal? }

  \item{weightESy}{ weight to use for the trait as an additional eigengene; should be between 0 and 1 }

  \item{getQValues}{ logical: should q-values be calculated? }

}
\details{
This function should be considered experimental. It takes into account both the "standard" and the network
measures of gene importance for the trait.
}
\value{
  datout = data.frame(p.Weighted, q.Weighted, Cor.Weighted,
        Z.Weighted, p.Standard, q.Standard, Cor.Standard, Z.Standard)
  Data frame reporting the following quantities for each given gene:

  \item{p.Weighted }{weighted p-value of association with the trait}

  \item{q.Weighted }{q-value (local FDR) calculated from \code{p.Weighted}}

  \item{cor.Weighted}{correlation of trait with gene expression weighted by a network term}

  \item{Z.Weighted}{ Fisher Z score of the weighted correlation}

  \item{p.Standard}{ standard Student p-value of association of the gene with the trait}

  \item{q.Standard}{ q-value (local FDR) calculated from \code{p.Standard}}

  \item{cor.Standard}{ correlation of gene with the trait}

  \item{Z.Standard}{ Fisher Z score of the standard correlation}

}
\author{ Steve Horvath }
\keyword{ misc}
