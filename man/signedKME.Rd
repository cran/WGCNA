\name{signedKME}
\alias{signedKME}
\title{ Signed eigengene-based connectivity }
\description{
  Calculation of (signed) eigengene-based connectivity, also known as module membership.
}
\usage{
signedKME(
  datExpr, datME, 
  outputColumnName = "kME", 
  corFnc = "cor", corOptions = "use = 'p'")
}
\arguments{
  \item{datExpr}{ a data frame containing the gene expression data. Rows correspond to samples and
columns to genes. Missing values are allowed and will be ignored. }

  \item{datME}{ a data frame containing module eigengenes. Rows  correspond to samples and columns to
module eigengenes. }

  \item{outputColumnName}{ a character string specifying the prefix of column names of the output. }

  \item{corFnc}{ character string specifying the function to be used to calculate co-expression
similarity. Defaults to Pearson correlation. Any function returning values between -1 and 1 can be used. }

  \item{corOptions}{ character string specifying additional arguments to be passed to the function given
by \code{corFnc}. Use \code{"use = 'p', method = 'spearman'"} to obtain Spearman correlation.   }

}
\details{
  Signed eigengene-based connectivity of a gene in a module is defined as the correlation of the gene
with the corresponding module eigengene.  The samples in \code{datExpr} and \code{datME} must be the
same.
}
\value{
  A data frame in which rows correspond to input genes and columns to module eigengenes, giving the
signed eigengene-based connectivity of each gene with respect to each eigengene.
}
\references{ 

Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24

Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression Network Analysis. PLoS Comput Biol
4(8): e1000117

}

\author{ Steve Horvath }

\keyword{ misc }% __ONLY ONE__ keyword per line
