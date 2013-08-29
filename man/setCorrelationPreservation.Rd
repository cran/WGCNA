\name{setCorrelationPreservation}
\alias{setCorrelationPreservation}
\title{ Summary correlation preservation measure }
\description{
  Given consensus eigengenes, the function calculates the average correlation preservation pair-wise for
all pairs of sets.
}
\usage{
setCorrelationPreservation(
   multiME, 
   setLabels, 
   excludeGrey = TRUE, greyLabel = "grey", 
   method = "absolute")
}
\arguments{

  \item{multiME}{ consensus module eigengenes in a multi-set format. A vector of lists with one list
corresponding to each set. Each list must contain a component \code{data} that is a data frame whose
columns are consensus module eigengenes. }

  \item{setLabels}{names to be used for the sets represented in \code{multiME}.}

  \item{excludeGrey}{logical: exclude the 'grey' eigengene from preservation measure?}

  \item{greyLabel}{module label corresponding to the 'grey' module. Usually this will be the
character string \code{"grey"} if the labels are colors, and the number 0 if the labels are numeric.}
  
  \item{method}{ character string giving the correlation preservation measure to use. Recognized values
are (unique abbreviations of) \code{"absolute"}, \code{"hyperbolic"}. }

}
\details{
  For each pair of sets, the function calculates the average preservation of correlation among the
eigengenes. Two preservation measures are available, the abosolute preservation (high if the two
correlations are similar and low if they are different), and the hyperbolically scaled preservation,
which de-emphasizes preservation of low correlation values.
}
\value{
  A data frame with each row and column corresponding to a set given in \code{multiME}, containing the
pairwise average correlation preservation values. Names and rownames are set to entries of \code{setLabels}. 
}
\references{ Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between
co-expression modules. BMC Systems Biology 2007, 1:54 }

\author{ Peter Langfelder }
\seealso{ 

\code{\link{multiSetMEs}} for module eigengene calculation;

\code{\link{plotEigengeneNetworks}} for eigengene network visualization.

}
\keyword{ misc }% __ONLY ONE__ keyword per line
