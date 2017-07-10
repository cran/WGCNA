\name{verboseBarplot}
\alias{verboseBarplot}
\title{ Barplot with error bars, annotated by Kruskal-Wallis or ANOVA p-value}
\description{
  Produce a barplot with error bars, annotated by Kruskal-Wallis or ANOVA p-value.
}
\usage{
verboseBarplot(x, g, 
               main = "", xlab = NA, ylab = NA, 
               cex = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
               color = "grey", numberStandardErrors = 1, 
               KruskalTest = TRUE, AnovaTest = FALSE, two.sided = TRUE, 
               addCellCounts=FALSE, horiz = FALSE, ylim = NULL, ...,
               addScatterplot = FALSE,
               pt.cex = 0.8, pch = 21, pt.col = "blue", pt.bg = "skyblue",
               randomSeed = 31425, jitter = 0.6,
               adjustYLim = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ numerical or binary vector of data whose group means are to be plotted }
  \item{g}{ a factor or a an object coercible to a factor giving the groups whose means are to be
calculated. }
  \item{main}{ main title for the plot.}
  \item{xlab}{ label for the x-axis. }
  \item{ylab}{ label for the y-axis. }
  \item{cex}{ character expansion factor for plot annotations. }
  \item{cex.axis}{ character expansion factor for axis annotations. }
  \item{cex.lab}{ character expansion factor for axis labels. }
  \item{cex.main}{ character expansion factor for the main title. }
  \item{color}{ a vector giving the colors of the bars in the barplot. }
  \item{numberStandardErrors}{ size of the error bars in terms of standard errors. See details. }
  \item{KruskalTest}{logical: should Kruskal-Wallis test be performed? See details. }
  \item{AnovaTest}{ logical: should ANOVA be performed? See details.  }
  \item{two.sided}{ logical: should the printed p-value be two-sided? See details. }
  \item{addCellCounts}{ logical: should counts be printed above each bar? }
  \item{horiz}{ logical: should the bars be drawn horizontally? }
  \item{ylim}{optional specification of the limits for the y axis. If not given, they will be determined
automatically.}
  \item{\dots}{ other parameters to function \code{\link{barplot}}.  }
  \item{addScatterplot}{logical: should a scatterplot of the data be overlaid? }
  \item{pt.cex}{character expansion factor for the points.}
  \item{pch}{shape code for the points.}
  \item{pt.col}{color for the points.}
  \item{pt.bg}{background color for the points.}
  \item{randomSeed}{integer random seed to make plots reproducible.}
  \item{jitter}{amount of random jitter to add to the position of the points along the x axis.}
  \item{adjustYLim}{logical: should the limits of the y axis be set so as to accomodate the individual points?
The adjustment is only carried out if input \code{ylim} is \code{NULL} and \code{addScatterplot} is
\code{TRUE}. In particular, if the user supplies \code{ylim}, it is not touched.}
}
\details{
  This function creates a barplot of a numeric variable (input \code{x}) across the levels of a grouping
variable (input \code{g}). The height of the bars equals the mean value of \code{x} across the
observations with a given level of \code{g}. By default, the barplot also shows plus/minus one standard
error. If you want only plus one standard error (not minus) choose \code{two.sided=TRUE}.  But the number
of standard errors can be determined with the input \code{numberStandardErrors}. For example, if you want
a 95\% confidence interval around the mean, choose \code{numberStandardErrors=2}. If you don't want any
standard errors set \code{numberStandardErrors=-1}.  The function also outputs the p-value of a Kruskal
Wallis test (Fisher test for binary input data), 
which is a non-parametric multi group comparison test. Alternatively, one can use Analysis
of Variance (Anova) to compute a p-value by setting \code{AnovaTest=TRUE}.  Anova is a generalization of
the Student t-test to multiple groups. In case of two groups, the Anova p-value equals the Student t-test
p-value. Anova should only be used if \code{x} follows a normal distribution. Anova also assumes
homoscedasticity (equal variances). The Kruskal Wallis test is often advantageous since it makes no
distributional assumptions.  Since the Kruskal Wallis test is based on the ranks of \code{x}, it is more
robust with regard to outliers. All p-values are two-sided. 
}
\value{
  None.
}
\author{ Steve Horvath, with contributions from Zhijin (Jean) Wu and Peter Langfelder}
\seealso{ \code{\link{barplot}} }
\examples{

   group=sample(c(1,2),100,replace=TRUE)

   height=rnorm(100,mean=group)

   par(mfrow=c(2,2))
   verboseBarplot(height,group, main="1 SE, Kruskal Test")

   verboseBarplot(height,group,numberStandardErrors=2, 
                  main="2 SE, Kruskal Test")

   verboseBarplot(height,group,numberStandardErrors=2,AnovaTest=TRUE, 
                  main="2 SE, Anova")

   verboseBarplot(height,group,numberStandardErrors=2,AnovaTest=TRUE, 
                  main="2 SE, Anova, only plus SE", two.sided=FALSE)

}
\keyword{ misc }
