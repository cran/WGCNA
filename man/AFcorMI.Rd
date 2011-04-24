\name{AFcorMI}
\alias{AFcorMI}
\title{Prediction of Weighted Mutual Information Adjacency Matrix by Correlation
}
\description{
AFcorMI computes a predicted weighted mutual information adjacency matrix from a given
correlation matrix.}
\usage{
AFcorMI(r, m)
}
\arguments{
  \item{r}{
a symmetric correlation matrix with values from -1 to 1.
}
  \item{m}{
number of observations from which the correlation was calcuated.
}
}
\details{
This function is a one-to-one prediction when we consider correlation as unsigned. The prediction
corresponds to the \code{AdjacencyUniversalVersion2} discussed in the help file for  the function
\code{\link{mutualInfoAdjacency}}. For more information
about the generation and features of the predicted mutual information adjacency, please refer to the function
\code{\link{mutualInfoAdjacency}}.
}
\value{
A matrix with the same size as the input correlation matrix, containing the predicted mutual information of
type  \code{AdjacencyUniversalVersion2}.
}
\author{
Steve Horvath, Lin Song, Peter Langfelder
}

\seealso{
\code{\link{mutualInfoAdjacency}}
}
\examples{
#Simulate a data frame datE which contains 5 columns and 50 observations
m=50
x1=rnorm(m)
r=.5; x2=r*x1+sqrt(1-r^2)*rnorm(m)
r=.3; x3=r*(x1-.5)^2+sqrt(1-r^2)*rnorm(m)
x4=rnorm(m)
r=.3; x5=r*x4+sqrt(1-r^2)*rnorm(m)
datE=data.frame(x1,x2,x3,x4,x5)
#calculate predicted AUV2
cor.data=cor(datE, use="p")
AUV2=AFcorMI(r=cor.data, m=nrow(datE))

}
\keyword{ misc }
