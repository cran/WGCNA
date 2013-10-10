\name{mutualInfoAdjacency}
\Rdversion{1.1}
\alias{mutualInfoAdjacency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate weighted adjacency matrices based on mutual information}
\description{
The function calculates different types of weighted adjacency matrices based on the mutual information
between vectors (corresponding to the columns of the input data frame datE).  The mutual information between
pairs of vectors is divided by an upper bound so that the resulting normalized measure lies between 0 and 1.
} 
\usage{
mutualInfoAdjacency(
   datE, 
   discretizeColumns = TRUE, 
   entropyEstimationMethod = "MM", 
   numberBins = NULL)
}
\arguments{
  \item{datE}{ \code{datE} is a data frame or matrix whose columns correspond to variables and whose rows correspond to measurements. For example, the columns may correspond to genes while the rows correspond to microarrays. The number of nodes in the mutual information network equals the number of columns of \code{datE}.
}
  \item{discretizeColumns}{ is a logical variable. If it is set to TRUE then the columns of \code{datE}  will be discretized into a user-defined number of bins (see \code{numberBins}).  
}
  \item{entropyEstimationMethod}{ takes a text string for specifying the entropy and mutual information estimation method. If \code{entropyEstimationMethod="MM"} then the Miller-Madow asymptotic bias corrected empirical estimator is used. 
 If \code{entropyEstimationMethod="ML"} the maximum likelihood estimator (also known as plug-in or empirical estimator) is used.
If  \code{entropyEstimationMethod="shrink"}, the shrinkage estimator of a Dirichlet probability distribution is used.
If  \code{entropyEstimationMethod="SG"}, the Schurmann-Grassberger estimator of the entropy of a Dirichlet probability distribution is used.
}
  \item{numberBins}{ is an integer larger than 0 which specifies how many bins are used for the discretization step. This argument is only relevant if \code{discretizeColumns} has been set to TRUE. By default \code{numberBins} is set to sqrt(m) where m is the number of samples, i.e. the number of rows of \code{datE}. Thus the default is \code{numberBins}=sqrt(nrow(datE)).   
}
}
\details{
The function inputs a data frame \code{datE} and outputs a list whose components correspond to different weighted network adjacency measures defined beteween the columns of \code{datE}. Make sure to install the following R packages \code{entropy}, \code{minet}, \code{infotheo} since 
the function \code{mutualInfoAdjacency} makes use of the \code{entropy} function from the R package \code{entropy} (Hausser and Strimmer 2008) and functions from the \code{minet} and \code{infotheo} package (Meyer et al 2008). 
A weighted network adjacency matrix is a symmetric matrix whose entries take on values between 0 and 1. Each weighted adjacency matrix contains scaled versions of the mutual information between the columns of the input data frame \code{datE}.
We assume that datE contains numeric values which will be discretized unless the user chooses the option \code{discretizeColumns=FALSE}.
The raw (unscaled) mutual information and entropy measures have units "nat", i.e. natural logarithms are used in their definition (base e=2.71..).  
Several mutual information estimation methods have been proposed in the literature (reviewed in Hausser and Strimmer 2008, Meyer et al 2008). 
While mutual information networks allows one to detect non-linear relationships between the columns of \code{datE}, they may overfit the data if relatively few observations are available. Thus, if the number of rows of \code{datE} is smaller than say 200, it may be better to fit a correlation using the function \code{adjacency}.
}
\value{
The function outputs a list with the following components:

 \item{ Entropy}{   is a vector whose components report entropy estimates of each column of \code{datE}. The natural logarithm (base e) is used in the definition. Using the notation from the Wikipedia entry (http://en.wikipedia.org/wiki/Mutual_information), this vector contains the values Hx where x corresponds to a column in \code{datE}.
 }
\item{MutualInformation}{   is a symmetric matrix whose entries contain the pairwise mutual information
measures between the columns of \code{datE}. The diagonal of the matrix \code{MutualInformation} equals
\code{Entropy}. In general, the entries of this matrix can be larger than 1, i.e. this is not an adjacency
matrix. Using the notation from the Wikipedia entry, this matrix contains the mutual information estimates
I(X;Y)  } 

\item{AdjacencySymmetricUncertainty}{   is a weighted adjacency matrix whose entries are based on the mutual
information. Using the notation from the Wikipedia entry, this matrix contains the mutual information
estimates \code{AdjacencySymmetricUncertainty}=2*I(X;Y)/(H(X)+H(Y)). Since I(X;X)=H(X), the diagonal
elements of \code{AdjacencySymmetricUncertainty} equal 1. In general the entries of this symmetric matrix
\code{AdjacencySymmetricUncertainty} lie between 0 and 1. 
} 

\item{AdjacencyUniversalVersion1}{ is a weighted adjacency matrix that is a simple function of the
\code{AdjacencySymmetricUncertainty}. Specifically, \code{AdjacencyUniversalVersion1=
AdjacencySymmetricUncertainty/(2- AdjacencySymmetricUncertainty)}. Note that f(x)= x/(2-x) is a
monotonically increasing function on the unit interval [0,1] whose values lie between 0 and 1. The reason
why we call it the universal adjacency is that dissUA=1-\code{AdjacencyUniversalVersion1} turns out to be
a universal distance function, i.e. it satisfies the properties of a distance (including the triangle
inequality) and it takes on a small value if any other distance measure takes on a small value (Kraskov et
al 2003). 
}

\item{AdjacencyUniversalVersion2}{ is a weighted adjacency matrix for which dissUAversion2=1-\code{AdjacencyUniversalVersion2} is also a universal distance measure. Using the notation from Wikipedia, the entries of the symmetric matrix AdjacencyUniversalVersion2 are defined as follows
\code{AdjacencyUniversalVersion2}=I(X;Y)/max(H(X),H(Y)).

}


}
\references{
	Hausser J, Strimmer K (2008) Entropy inference and the James-Stein 	estimator, with application to nonlinear gene association networks. See 	http://arxiv.org/abs/0811.3579

Patrick E. Meyer, Frederic Lafitte, and Gianluca Bontempi. minet: A R/Bioconductor Package for Inferring Large Transcriptional Networks Using Mutual Information. BMC Bioinformatics, Vol 9, 2008

Kraskov A, Stoegbauer H, Andrzejak RG, Grassberger P (2003) Hierarchical Clustering Based on Mutual Information. ArXiv q-bio/0311039

}
\author{
Steve Horvath, Lin Song, Peter Langfelder
}

\seealso{
\code{\link{adjacency}}
}
\examples{

# Load requisite packages. These packages are considered "optional", 
# so WGCNA does not load them automatically.

if (require(infotheo, quietly = TRUE) && 
    require(minet, quietly = TRUE) && 
    require(entropy, quietly = TRUE))
{
  # Example can be executed.
  #Simulate a data frame datE which contains 5 columns and 50 observations
  m=50
  x1=rnorm(m)
  r=.5; x2=r*x1+sqrt(1-r^2)*rnorm(m)
  r=.3; x3=r*(x1-.5)^2+sqrt(1-r^2)*rnorm(m)
  x4=rnorm(m)
  r=.3; x5=r*x4+sqrt(1-r^2)*rnorm(m)
  datE=data.frame(x1,x2,x3,x4,x5)
  
  #calculate entropy, mutual information matrix and weighted adjacency 
  # matrices based on mutual information.
  MIadj=mutualInfoAdjacency(datE=datE)
} else 
  printFlush(paste("Please install packages infotheo, minet and entropy",
                   "before running this example."));

}

\keyword{ misc }
