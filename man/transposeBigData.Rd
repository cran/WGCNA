\name{transposeBigData}
\alias{transposeBigData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transpose a big matrix or data frame
}
\description{
This transpose command partitions a big matrix (or data frame) into blocks and applies the t() function to
each block separately.  
}
\usage{
transposeBigData(x, blocksize = 20000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix or data frame
}
  \item{blocksize}{a positive integer larger than 1, which determines the block size. Default is 20k.
}
}
\details{
Assume you have a very large matrix with say 500k columns. In this case, the standard transpose function of
R \code{t()} can take a long time. Solution: Split the original matrix into sub-matrices by dividing the
columns into blocks. Next apply  \code{t()} to each sub-matrix. The same holds if the large matrix contains
a large number of rows.  The function \code{transposeBigData} automatically checks whether the large matrix
contains more rows or more columns. If the number of columns is larger than or equal to the number of rows
then the block wise splitting will be applied to columns otherwise to the rows.  } 

\value{A matrix or data frame (depending on the input \code{x} ) which is the transpose of \code{x}. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Any linear algebra book will explain the transpose.
}
\author{
Steve Horvath, UCLA
}
\note{
This function can be considered a wrapper of \code{\link{t}()}
}

\seealso{
The standard function \code{\link{t}} .
}
\examples{
x=data.frame(matrix(1:10000,nrow=4,ncol=2500))
dimnames(x)[[2]]=paste("Y",1:2500,sep="")
xTranspose=transposeBigData(x)
x[1:4,1:4]
xTranspose[1:4,1:4]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
