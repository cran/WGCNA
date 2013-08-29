\name{populationMeansInAdmixture}
\alias{populationMeansInAdmixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the population-specific mean values in an admixed population.
}
\description{Uses the expression values from an admixed population and estimates of the proportions of
sub-populations to estimate the population specific mean values. For example, this function can be used to
estimate the cell type specific mean gene expression values based on expression values from a mixture of
cells. The method is described in Shen-Orr et al (2010) where it was used to estimate cell type specific
gene expression levels based on a mixture sample. }
\usage{
 populationMeansInAdmixture(
    datProportions, datE.Admixture, 
    scaleProportionsTo1 = TRUE,
    scaleProportionsInCelltype = TRUE, 
    setMissingProportionsToZero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datProportions}{a matrix of non-negative numbers (ideally proportions) where the rows correspond to the samples (rows of \code{datE.Admixture}) and the columns correspond to the sub-populations of the mixture. The function calculates a mean expression value for each column of \code{datProportions}.
  Negative entries in datProportions lead to an error message. But the rows of
  \code{datProportions} do not have to sum to 1, see the argument \code{scaleProportionsTo1}.
}
  \item{datE.Admixture}{a matrix of numbers. The rows correspond to samples (mixtures of populations). The columns contain the variables (e.g. genes) for which the means should be estimated.  
}
  \item{scaleProportionsTo1}{logical. If set to TRUE (default) then the proportions in each row of \code{datProportions} are scaled so that they sum to 1, i.e. datProportions[i,]=datProportions[i,]/max(datProportions[i,]).
  In general, we recommend to set it to TRUE.
}
 \item{scaleProportionsInCelltype}{
 logical. If set to TRUE (default) then the proportions in each cell types are recaled and make the mean to
 0.
}
  \item{setMissingProportionsToZero}{logical. Default is FALSE. If set to TRUE then it sets missing values in \code{datProportions} to zero.
}
}
\details{The function outputs a matrix of coefficients resulting from fitting a regression model. If the proportions sum to 1, then i-th row of the output matrix reports the coefficients  of the following model
\code{lm(datE.Admixture[,i]~.-1,data=datProportions)}. Aside, the minus 1 in the formula indicates that no intercept term will be fit.
Under certain assumptions, the coefficients can be interpreted as the mean expression values in the sub-populations (Shen-Orr  2010).
}
\value{ a numeric matrix whose rows correspond to the columns of \code{datE.Admixture} (e.g. to genes) and whose columns correspond to the columns of \code{datProportions} (e.g. sub populations or cell types).
}
\references{
Shen-Orr SS, Tibshirani R, Khatri P, Bodian DL, Staedtler F, Perry NM,
Hastie T, Sarwal MM, Davis MM, Butte AJ (2010)  
Cell type-specific gene expression differences in complex
tissues.  Nature Methods,  vol 7 no.4 
}
\author{
Steve Horvath, Chaochao Cai
}
\note{
This can be considered a wrapper of the \code{lm} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(1)
# this is the number of complex (mixed) tissue samples, e.g. arrays
m=10
# true count data (e.g. pure cells in the mixed sample)
datTrueCounts=as.matrix(data.frame(TrueCount1=rpois(m,lambda=16),
TrueCount2=rpois(m,lambda=8),TrueCount3=rpois(m,lambda=4),
TrueCount4=rpois(m,lambda=2)))
no.pure=dim(datTrueCounts)[[2]]

# now we transform the counts into proportions
divideBySum=function(x) t(x)/sum(x)
datProportions= t(apply(datTrueCounts,1,divideBySum))
dimnames(datProportions)[[2]]=paste("TrueProp",1:dim(datTrueCounts)[[2]],sep=".")

# number of genes that are highly expressed in each pure population
no.genesPerPure=rep(5, no.pure)
no.genes= sum(no.genesPerPure)
GeneIndicator=rep(1:no.pure, no.genesPerPure)
# true mean values of the genes in the pure populations
# in the end we hope to estimate them from the mixed samples
datTrueMeans0=matrix( rnorm(no.genes*no.pure,sd=.3), nrow= no.genes,ncol=no.pure)
for (i in 1:no.pure ){
datTrueMeans0[GeneIndicator==i,i]= datTrueMeans0[GeneIndicator==i,i]+1
}
dimnames(datTrueMeans0)[[1]]=paste("Gene",1:dim(datTrueMeans0)[[1]],sep="." )
dimnames(datTrueMeans0)[[2]]=paste("MeanPureCellType",1:dim(datTrueMeans0)[[2]],
                                   sep=".")
# plot.mat(datTrueMeans0)
# simulate the (expression) values of the admixed population samples

noise=matrix(rnorm(m*no.genes,sd=.1),nrow=m,ncol= no.genes)
datE.Admixture= as.matrix(datProportions) \%*\% t(datTrueMeans0) + noise
dimnames(datE.Admixture)[[1]]=paste("MixedTissue",1:m,sep=".")

datPredictedMeans=populationMeansInAdmixture(datProportions,datE.Admixture)

par(mfrow=c(2,2))
for (i in 1:4 ){
verboseScatterplot(datPredictedMeans[,i],datTrueMeans0[,i],
xlab="predicted mean",ylab="true mean",main="all populations")
abline(0,1)
}

#assume we only study 2 populations (ie we ignore the others)
selectPopulations=c(1,2)
datPredictedMeansTooFew=populationMeansInAdmixture(datProportions[,selectPopulations],
                                                   datE.Admixture)

par(mfrow=c(2,2))
for (i in 1:length(selectPopulations) ){
verboseScatterplot(datPredictedMeansTooFew[,i],datTrueMeans0[,i],
xlab="predicted mean",ylab="true mean",main="too few populations")
abline(0,1)
}

#assume we erroneously add a population
datProportionsTooMany=data.frame(datProportions,WrongProp=sample(datProportions[,1]))
datPredictedMeansTooMany=populationMeansInAdmixture(datProportionsTooMany,
                                 datE.Admixture)

par(mfrow=c(2,2))
for (i in 1:4 ){
  verboseScatterplot(datPredictedMeansTooMany[,i],datTrueMeans0[,i],
  xlab="predicted mean",ylab="true mean",main="too many populations")
  abline(0,1)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{misc}% __ONLY ONE__ keyword per line
