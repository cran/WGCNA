\name{plotModuleSignificance}
\alias{plotModuleSignificance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Barplot of module significance }
\description{
  Plot a barplot of gene significance. 
}
\usage{
plotModuleSignificance(
  geneSignificance, 
  colors, 
  boxplot = FALSE, 
  main = "Gene significance across modules,", 
  ylab = "Gene Significance", ...)
}
\arguments{
  \item{geneSignificance}{ a numeric vector giving gene significances. }

  \item{colors}{ a character vector specifying module assignment for the genes whose significance is
given in \code{geneSignificance }. The modules should be labeled by colors. }

  \item{boxplot}{ logical: should a boxplot be produced instead of a barplot? }

  \item{main}{ main title for the plot. }
  \item{ylab}{ y axis label for the plot. }
  \item{\dots}{ other graphical parameters to \code{\link{plot}}. }
}
\details{
  Given individual gene significances and their module assigment, the function calculates the module
significance for each module as the average gene significance of the genes within the module. The result
is plotted in a barplot or boxplot form. Each bar or box is labeled by the corresponding module color.
}
\value{
  None.
}
\references{ 

Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression
Network Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17

Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24

}

\author{ Steve Horvath }

\seealso{ \code{\link{barplot}}, \code{\link{boxplot}} }
\keyword{ hplot }% __ONLY ONE__ keyword per line
\keyword{ misc }% __ONLY ONE__ keyword per line
