\name{hubGeneSignificance}
\alias{hubGeneSignificance}
\title{ Hubgene significance }
\description{
  Calculate approximate hub gene significance for all modules in network.
}
\usage{
hubGeneSignificance(datKME, GS)
}
\arguments{
  \item{datKME}{ a data frame (or a matrix-like object) containing eigengene-based connectivities of all
genes in the network. }
  \item{GS}{ a vector with one entry for every gene containing its gene significance. }
}
\details{
  In \code{datKME} rows correspond to genes and columns to modules. 
}
\value{
  A vector whose entries are the hub gene significances for each module.
}
\references{ Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology
2007, 1:24 }
\author{ Steve Horvath }
\keyword{ misc }% __ONLY ONE__ keyword per line
