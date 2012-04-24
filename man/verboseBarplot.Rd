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
               horiz = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ numerical vector of data whose group means are to be plotted }
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
  \item{horiz}{ logical: should the bars be drawn horizontally? }
  \item{\dots}{ other parameters to function \code{\link{barplot}}  }
}
\details{
  This function creates a barplot of a numeric variable (input \code{x}) across the levels of a grouping
variable (input \code{g}). The height of the bars equals the mean value of \code{x} across the
observations with a given level of \code{g}. By default, the barplot also shows plus/minus one standard
error. If you want only plus one standard error (not minus) choose \code{two.sided=TRUE}.  But the number
of standard errors can be determined with the input \code{numberStandardErrors}. For example, if you want
a 95\% confidence interval around the mean, choose \code{numberStandardErrors=2}. If you don't want any
standard errors set \code{numberStandardErrors=-1}.  The function also outputs the p-value of a Kruskal
Wallis test, which is a non-parametric multi group comparison test. Alternatively, one can use Analysis
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
\author{ Steve Horvath }
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
